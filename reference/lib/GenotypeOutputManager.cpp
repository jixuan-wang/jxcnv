#include "GenotypeOutputManager.hpp"
#include "SampleGenotypeQualsCalculator.hpp"
#include "CNVmodelParams.hpp"
#include "AlleleQuals.hpp"
#include "CopyNumberVariant.hpp"
#include "AuxiliaryCNVoutputManager.hpp"
#include "xhmmOutputManager.hpp"

#include <params.hpp>
#include <Model.hpp>
#include <Data.hpp>
#include <Timer.hpp>
#include <utils.hpp>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
using namespace std;

#define USE_DB_NO_STDIN "TMP_GT_DB"
#define DB_SUFFIX ".db"
#define DB_INDEX_SUFFIX ".dbi"

const string XHMM::GenotypeOutputManager::MISSING_VAL = ".";

XHMM::GenotypeOutputManager::GenotypeOutputManager(string outFile, const DiscoverOutputManager::IntervalToThreshMap* intervals, AuxiliaryCNVoutputManager* auxOutput, const BaseReal genotypeQualThresholdWhenNoExact, const bool genotypeSubsegments, const uint maxTargetsInSubsegment, ReadDepthMatrixLoader* origDataLoader, BaseReal maxScoreVal, int scoreDecimalPrecision)
: CNVoutputManager(origDataLoader),
  _outStream(HMM_PP::utils::getOstreamWriterFromFile(outFile)),
  _intervals(new IntervalToPrecisionThreshMap()),
  _auxOutput(auxOutput),
  _genotypeQualThresholdWhenNoExact(genotypeQualThresholdWhenNoExact),
  _genotypeSubsegments(genotypeSubsegments), _maxTargetsInSubsegment(maxTargetsInSubsegment),
  _maxScoreVal(maxScoreVal), _scoreDecimalPrecision(scoreDecimalPrecision),
  _samples(new vector<string>()), _allGenotypes(NULL),
  _vcfTransposerFile((outFile != HMM_PP::utils::STD_STREAM ? outFile : USE_DB_NO_STDIN)),
  _vcfTransposer(NULL) {

	if (_maxScoreVal < _genotypeQualThresholdWhenNoExact)
		_maxScoreVal = _genotypeQualThresholdWhenNoExact;

	for (DiscoverOutputManager::IntervalToThreshMap::const_iterator it = intervals->begin(); it != intervals->end(); ++it) {
		BaseReal useThresh = _genotypeQualThresholdWhenNoExact;
		if (it->second.first != NULL) {
			useThresh = *(it->second.first);
			delete it->second.first;
		}
		(*_intervals)[it->first] = ThreshAndSamples(Genotype::RealThresh(_scoreDecimalPrecision, useThresh), it->second.second);
#if defined(DEBUG)
		cerr << "Will genotype CNV " << it->first << " with a quality of " << (*_intervals)[it->first].first << " (" << useThresh << ")" << endl;
#endif
	}
	delete intervals;

	if (_outStream == NULL)
		throw new Exception("Cannot give empty file name for output");

	initOstream((*_outStream)());
}

XHMM::GenotypeOutputManager::~GenotypeOutputManager() {
	if (_vcfTransposer != NULL) {
		printVCFgenotypesAllIntervals();
		delete _vcfTransposer;
	}

	delete _outStream;
	if (_intervals != NULL)
		delete _intervals;
	delete _samples;

	for (map<Interval, CopyNumberVariant*>::const_iterator it = _allGenotypes->begin(); it != _allGenotypes->end(); ++it)
		delete it->second;
	delete _allGenotypes;
}

void XHMM::GenotypeOutputManager::setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader) {
	CNVoutputManager::setModelAndReadDepthLoader(model, dataLoader);
	_allGenotypes = new map<Interval, CopyNumberVariant*>();

	bool storeGenotypes = false;

	// TODO: This all-by-all intervals-by-_dataLoader search can perhaps be done more efficiently using a variation of interval trees:
	for (IntervalToPrecisionThreshMap::const_iterator intervIt = _intervals->begin(); intervIt != _intervals->end(); ++intervIt) {
		const Interval& interv = intervIt->first;

		HMM_PP::ullintPair* rangeForGenotyping = NULL;
		if (_dataLoader->hasTargetsForChr(interv.getChr())) {
			for (uint t = _dataLoader->getChrStartTargetIndex(interv.getChr()); t <= _dataLoader->getChrStopTargetIndex(interv.getChr()); ++t) {
				const Interval& targT = _dataLoader->getTarget(t);

				if (interv.strictlyOverlaps(targT)) {
					if (rangeForGenotyping == NULL)
						rangeForGenotyping = new HMM_PP::ullintPair(t, t);

					// extend the range to t:
					rangeForGenotyping->second = t;
				}
				else if (rangeForGenotyping != NULL)
					break; // already found a valid interval, yet targT does not overlap, so nothing else will [since _dataLoader has sorted targets]
			}
		}

		if (rangeForGenotyping == NULL)
			XHMM::xhmmOutputManager::printWarning("Cannot genotype CNV in region " + interv.intervalString() + " since it has no overlapping targets - setting all data values to missing");

		CopyNumberVariant* cnv = new CopyNumberVariant(CNVmodelParams::NON_DIPLOID_TYPES, rangeForGenotyping, _dataLoader, intervIt->second.first, storeGenotypes, intervIt->second.second);
		delete rangeForGenotyping;
		(*_allGenotypes)[interv] = cnv;
	}

	delete _intervals;
	_intervals = NULL;

	if (_genotypeSubsegments) {
		// Add all sub-segment of _allGenotypes overlapping _maxTargetsInSubsegment or fewer targets:
		map<Interval, CopyNumberVariant*>* addGenotypes = new map<Interval, CopyNumberVariant*>();
		for (map<Interval, CopyNumberVariant*>::const_iterator intervIt = _allGenotypes->begin(); intervIt != _allGenotypes->end(); ++intervIt) {
			const CopyNumberVariant* cnv = intervIt->second;
			const TargetRange* tRange = cnv->getTargetRange();

			if (tRange != NULL) {
				for (uint t1 = tRange->getStartTargIndex(); t1 <= tRange->getStopTargIndex(); ++t1) {
					for (uint t2 = t1; t2 <= min(tRange->getStopTargIndex(), t1 + _maxTargetsInSubsegment - 1); ++t2) {
						const Interval subSegment = _dataLoader->getTarget(t1) + _dataLoader->getTarget(t2);
						if (_allGenotypes->find(subSegment) != _allGenotypes->end())
							continue; // don't add a segment given in the input

						BaseReal useThresh = max(_genotypeQualThresholdWhenNoExact, cnv->getGenotypingThreshold().getThreshold());
						map<Interval, CopyNumberVariant*>::iterator findAddIt = addGenotypes->find(subSegment);
						if (findAddIt != addGenotypes->end()) {
							CopyNumberVariant* addCNV = findAddIt->second;
							if (useThresh < addCNV->getGenotypingThreshold().getThreshold()) { // want to be as permissive as was originally possible
								delete addCNV;
								findAddIt = addGenotypes->end();
							}
						}

						if (findAddIt == addGenotypes->end()) {
							HMM_PP::ullintPair subRange(t1, t2);
							(*addGenotypes)[subSegment] = new CopyNumberVariant(CNVmodelParams::NON_DIPLOID_TYPES, &subRange, _dataLoader, Genotype::RealThresh(_scoreDecimalPrecision, useThresh), storeGenotypes, NULL);
						}
					}
				}
			}
		}

		_allGenotypes->insert(addGenotypes->begin(), addGenotypes->end());
		delete addGenotypes;
		cerr << "Will genotype a total of " << _allGenotypes->size() << " unique intervals and sub-intervals" << endl;
	}

	printVCFheader();

	_vcfTransposer = new OnDiskMatrixTransposer(_vcfTransposerFile, _allGenotypes->size());
}

void XHMM::GenotypeOutputManager::printDataType(const HMM_PP::InferredDataFit* idf) {
	HMM_PP::Data* origData = CNVoutputManager::getOrigDataForNextDataFitOutput(idf);

	SampleGenotypeQualsCalculator* genotypeCalc = new SampleGenotypeQualsCalculator(idf, _model, _maxScoreVal, _dataLoader);

	string samp = genotypeCalc->getSample();
	_samples->push_back(samp);

	list<string>* sampleGenotypes = new list<string>();

	HMM_PP::Timer* t = NULL;
	if (idf->isCalcTimes())
		t = new HMM_PP::Timer();

	for (map<Interval, CopyNumberVariant*>::iterator intervIt = _allGenotypes->begin(); intervIt != _allGenotypes->end(); ++intervIt) {
		CopyNumberVariant* cnv = intervIt->second;
		const TargetRange* rng = cnv->getTargetRange();

		const HMM_PP::ModelParams::DataVal* meanOrigRD = NULL;
		if (origData != NULL && rng != NULL)
			meanOrigRD = _model->calcRepresentativeDataVal(origData, rng->getStartTargIndex(), rng->getStopTargIndex());

		const Genotype::RealThresh& thresh = cnv->getGenotypingThreshold();
		Genotype* gt = new Genotype(genotypeCalc, rng, CNVmodelParams::DIPLOID, CNVmodelParams::NON_DIPLOID_TYPES, meanOrigRD, &thresh);
		sampleGenotypes->push_back(generateGenotypeString(gt, cnv));
		cnv->addGenotype(samp, gt);

		if (_auxOutput != NULL && rng != NULL)
			_auxOutput->printCNVtargets(genotypeCalc, rng->getStartTargIndex(), rng->getStopTargIndex(), CNVmodelParams::ALL_CN_TYPES, origData);
	}

	if (idf->isCalcTimes()) {
		idf->addCalcTime("GT", t->getDuration());
		delete t;
	}

	if (origData != NULL)
		delete origData;
	delete genotypeCalc;

	if (idf->isCalcTimes())
		t = new HMM_PP::Timer();
	_vcfTransposer->cacheNextRowToDBandDelete(sampleGenotypes);

	if (idf->isCalcTimes()) {
		idf->addCalcTime("DB", t->getDuration());
		delete t;
	}
}

const string XHMM::GenotypeOutputManager::HEADER_ROW_PREFIX = "##";

const string XHMM::GenotypeOutputManager::VCF_COLUMNS_ARRAY[] = {
		"#CHROM",
		"POS",
		"ID",
		"REF",
		"ALT",
		"QUAL",
		"FILTER",
		"INFO",
		"FORMAT"
};
list<string> XHMM::GenotypeOutputManager::VCF_COLUMNS = HMM_PP::arrayToList(XHMM::GenotypeOutputManager::VCF_COLUMNS_ARRAY);

void XHMM::GenotypeOutputManager::printVCFheader() {
	ostream& outStream = (*_outStream)();

	const uint nAlt = CNVmodelParams::NON_DIPLOID_TYPES.size();
	const uint numTotalAlleles = CNVmodelParams::ALL_CN_TYPES.size();

	outStream
	<< HEADER_ROW_PREFIX << "fileformat=VCFv4.1" << endl
	<< HEADER_ROW_PREFIX << "ALT=<ID=" << _model->getModelParams()->state(CNVmodelParams::DIPLOID) << ",Description=\"Diploid copy number\">" << endl
	<< HEADER_ROW_PREFIX << "ALT=<ID=CNV,Description=\"Copy Number Polymorphism\">" << endl
	<< HEADER_ROW_PREFIX << "ALT=<ID=" << _model->getModelParams()->state(CNVmodelParams::DEL) << ",Description=\"Deletion\">" << endl
	<< HEADER_ROW_PREFIX << "ALT=<ID=" << _model->getModelParams()->state(CNVmodelParams::DUP) << ",Description=\"Duplication\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=AC,Number=" << nAlt << ",Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=AF,Number=" << nAlt << ",Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=TPOS,Number=1,Type=Integer,Description=\"Start coordinate of target used to genotype this variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=TEND,Number=1,Type=Integer,Description=\"End coordinate of target used to genotype this variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=NUMT,Number=1,Type=Integer,Description=\"Number of targets used to genotype this variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=GQT,Number=1,Type=Float,Description=\"CNV-specific genotyping quality threshold, calculated as the minimal " << DiscoverOutputManager::EXACT_QUAL << " of discovered CNVs\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=PREVTARGSTART,Number=1,Type=Integer,Description=\"Start coordinate of target preceding the first target used to genotype this variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=PREVTARGEND,Number=1,Type=Integer,Description=\"End coordinate of target preceding the first target used to genotype this variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=POSTTARGSTART,Number=1,Type=Integer,Description=\"Start coordinate of target following the last target used to genotype this variant\">" << endl
	<< HEADER_ROW_PREFIX << "INFO=<ID=POSTTARGEND,Number=1,Type=Integer,Description=\"End coordinate of target following the last target used to genotype this variant\">" << endl;

	outStream
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=NDQ,Number=1,Type=Float,Description=\"Phred-scaled quality of =N=on =D=iploidy\">" << endl
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=DQ,Number=1,Type=Float,Description=\"Phred-scaled quality of =D=iploidy\">" << endl
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=EQ,Number=" << nAlt << ",Type=Float,Description=\"Phred-scaled qualities of =E=xact CNV event of allele types, in order given in ALT field\">" << endl
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=SQ,Number=" << nAlt << ",Type=Float,Description=\"Phred-scaled qualities of =S=ome CNV event of allele types, in order given in ALT field\">" << endl
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=NQ,Number=" << nAlt << ",Type=Float,Description=\"Phred-scaled qualities of =N=o CNV event of allele types, in order given in ALT field\">" << endl
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=LQ,Number=" << nAlt << ",Type=Float,Description=\"Phred-scaled qualities of =L=eft breakpoint of CNV event of allele types, in order given in ALT field\">" << endl
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=RQ,Number=" << nAlt << ",Type=Float,Description=\"Phred-scaled qualities of =R=ight breakpoint of CNV event of allele types, in order given in ALT field\">" << endl;

	outStream
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=PL,Number=" << numTotalAlleles << ",Type=Float,Description=\"Normalized, Phred-scaled relative likelihoods for ";
	for (uint i = 0; i < CNVmodelParams::ALL_CN_TYPES.size(); ++i) {
		if (i > 0)
			outStream << ",";
		outStream
		<< _model->getModelParams()->state(CNVmodelParams::ALL_CN_TYPES[i]);
	}
	outStream
	<< " genotypes, capped at " << setprecision(_scoreDecimalPrecision) << Genotype::MAX_PL << setprecision(DEFAULT_PRECISION) << "\">" << endl;

	outStream
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=RD,Number=1,Type=Float,Description=\"Mean Read Depth over region\">" << endl;

	if (hasOrigData())
		outStream
		<< HEADER_ROW_PREFIX << "FORMAT=<ID=ORD,Number=1,Type=Float,Description=\"Mean Original (unnormalized) Read Depth over region\">" << endl;

	outStream
	<< HEADER_ROW_PREFIX << "FORMAT=<ID=DSCVR,Number=1,Type=Character,Description=\"Was this CNV locus discovered in this sample? (Y or N)\">" << endl;
}

void XHMM::GenotypeOutputManager::printSampleHeaderLine(ostream& outStream, const vector<string>* samples) {
	// Print the SAMPLE header line:
	for (list<string>::const_iterator i = GenotypeOutputManager::VCF_COLUMNS.begin(); i != GenotypeOutputManager::VCF_COLUMNS.end(); ++i) {
		if (i != GenotypeOutputManager::VCF_COLUMNS.begin())
			outStream << '\t';
		outStream << *i;
	}

	for (vector<string>::const_iterator sampIt = samples->begin(); sampIt != samples->end(); ++sampIt)
		outStream << '\t' << *sampIt;
	outStream << endl;
}

void XHMM::GenotypeOutputManager::printVCFgenotypesAllIntervals() {
	ostream& outStream = (*_outStream)();

	printSampleHeaderLine(outStream, _samples);

	// Print each CNV:
	for (map<Interval, CopyNumberVariant*>::iterator intervIt = _allGenotypes->begin(); intervIt != _allGenotypes->end(); ++intervIt) {
		const Interval& interv = intervIt->first;
		const CopyNumberVariant* cnv = intervIt->second;

		const TargetRange* rng = cnv->getTargetRange();
		const Interval* regionUsed = cnv->getMergedTargets();
		const uint numTargsUsed = rng == NULL ? 0 : rng->getNumTargets();

		outStream
		<< interv.getChr() << '\t'
		<< interv.getBp1() << '\t' // TODO: Should we subtract 1 here, so that it's positioned one base BEFORE the event starts??? [I think that's what Bob Handsaker once said...]
		<< interv << '\t'
		<< "<" << _model->getModelParams()->state(CNVmodelParams::DIPLOID) << ">" << '\t'; // reference allele

		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			if (i > 0)
				outStream << ",";
			outStream
			<< "<" << _model->getModelParams()->state(CNVmodelParams::NON_DIPLOID_TYPES[i]) << ">";
		}
		outStream
		<< '\t';

		outStream
		<< MISSING_VAL << '\t'
		<< MISSING_VAL << '\t';

		uint an = cnv->getCalledSampleCount();
		vector<uint>* acVec = new vector<uint>(CNVmodelParams::NON_DIPLOID_TYPES.size());
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i)
			(*acVec)[i] = cnv->getNonRefTypeCount(CNVmodelParams::NON_DIPLOID_TYPES[i]);

		outStream
		<< "AC=";
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			if (i > 0)
				outStream << ",";
			outStream
			<< (*acVec)[i];
		}

		outStream
		<< ";AF=";
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			BaseReal af = 0;
			if (an > 0)
				af = (*acVec)[i] / (BaseReal) an;

			if (i > 0)
				outStream << ",";
			outStream
			<< af;
		}

		delete acVec;

		outStream
		<< ";AN=" << an;

		outStream
		<< ";END=" << interv.getBp2()
		<< ";IMPRECISE;SVLEN=" << interv.span()
		<< ";SVTYPE=" << "CNV";

		if (regionUsed != NULL)
			outStream
			<< ";TPOS=" << regionUsed->getBp1()
			<< ";TEND=" << regionUsed->getBp2();

		outStream
		<< ";NUMT=" << numTargsUsed
		<< ";GQT=" << cnv->getGenotypingThreshold();

		const Interval* prevTarget = cnv->getPreviousTarget();
		if (prevTarget != NULL && regionUsed != NULL && prevTarget->getChr() == regionUsed->getChr())
			outStream
			<< ";PREVTARGSTART=" << prevTarget->getBp1()
			<< ";PREVTARGEND=" << prevTarget->getBp2();

		const Interval* postTarget = cnv->getPostTarget();
		if (postTarget != NULL && regionUsed != NULL && postTarget->getChr() == regionUsed->getChr())
			outStream
			<< ";POSTTARGSTART=" << postTarget->getBp1()
			<< ";POSTTARGEND=" << postTarget->getBp2();

		outStream
		<< '\t';

		outStream
		<< "GT"
		<< ":" << "NDQ"
		<< ":" << "DQ"
		<< ":" << "EQ"
		<< ":" << "SQ"
		<< ":" << "NQ"
		<< ":" << "LQ"
		<< ":" << "RQ"
		<< ":" << "PL"
		<< ":" << "RD";

		if (hasOrigData())
			outStream
			<< ":" << "ORD";

		outStream
		<< ":" << "DSCVR";

		// Print the Genotype information for each sample at this CNV event:
		PullOutputDataToMap* mapData = new PullOutputDataToMap(_samples);
		_vcfTransposer->nextTransposedRowFromDB(mapData);
		SampleToGenotypeStrings* sampGts = mapData->getSampleGenotypes();

		for (vector<string>::const_iterator sampIt = _samples->begin(); sampIt != _samples->end(); ++sampIt) {
			const string& sample = *sampIt;

			outStream
			<< '\t'
			<< (*sampGts)[sample];
		}
		outStream << endl;

		delete mapData;
	}
}

string XHMM::GenotypeOutputManager::generateGenotypeString(const Genotype* gt, const CopyNumberVariant* cnv) const {
	stringstream gtStr;
	initOstream(gtStr);

	if (gt->isNoCall())
		gtStr << MISSING_VAL;
	else if (gt->isReference()) // the 0-th allele
		gtStr << 0;
	else // gt->isNonReference()
		gtStr << (cnv->getNonRefTypeIndex(gt->getNonRefCall()->getCNVtype()) + 1); // Add 1 for REF allele [whose index is 0]

	gtStr
	<< setprecision(_scoreDecimalPrecision);

	gtStr
	<< ":";
	const BaseReal* NDQ = gt->getNonDiploidScore();
	if (NDQ != NULL)
		gtStr
		<< *NDQ;
	else
		gtStr
		<< MISSING_VAL;

	gtStr
	<< ":";
	const BaseReal* DQ = gt->getDiploidScore();
	if (DQ != NULL)
		gtStr
		<< *DQ;
	else
		gtStr
		<< MISSING_VAL;

	gtStr
	<< ":";
	if (gt->getNumNonRefAlleles() == 0)
		gtStr
		<< MISSING_VAL;
	else {
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			if (i > 0)
				gtStr << ",";
			gtStr
			<< gt->getNonRefAllele(CNVmodelParams::NON_DIPLOID_TYPES[i])->getHaveExactEventScore();
		}
	}

	gtStr
	<< ":";
	if (gt->getNumNonRefAlleles() == 0)
		gtStr
		<< MISSING_VAL;
	else {
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			if (i > 0)
				gtStr << ",";
			gtStr
			<< gt->getNonRefAllele(CNVmodelParams::NON_DIPLOID_TYPES[i])->getHaveSomeEventScore();
		}
	}

	gtStr
	<< ":";
	if (gt->getNumNonRefAlleles() == 0)
		gtStr
		<< MISSING_VAL;
	else {
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			if (i > 0)
				gtStr << ",";
			gtStr
			<< gt->getNonRefAllele(CNVmodelParams::NON_DIPLOID_TYPES[i])->getDontHaveAnyEventScore();
		}
	}

	gtStr
	<< ":";
	if (gt->getNumNonRefAlleles() == 0)
		gtStr
		<< MISSING_VAL;
	else {
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			if (i > 0)
				gtStr << ",";
			gtStr
			<< gt->getNonRefAllele(CNVmodelParams::NON_DIPLOID_TYPES[i])->getStartStopScores().first;
		}
	}

	gtStr
	<< ":";
	if (gt->getNumNonRefAlleles() == 0)
		gtStr
		<< MISSING_VAL;
	else {
		for (uint i = 0; i < CNVmodelParams::NON_DIPLOID_TYPES.size(); ++i) {
			if (i > 0)
				gtStr << ",";
			gtStr
			<< gt->getNonRefAllele(CNVmodelParams::NON_DIPLOID_TYPES[i])->getStartStopScores().second;
		}
	}

	gtStr
	<< ":";
	if (gt->getNumNonRefAlleles() == 0)
		gtStr
		<< MISSING_VAL;
	else {
		for (uint i = 0; i < CNVmodelParams::ALL_CN_TYPES.size(); ++i) {
			if (i > 0)
				gtStr << ",";
			gtStr
			<< *(gt->getPL(CNVmodelParams::ALL_CN_TYPES[i]));
		}
	}

	gtStr
	<< setprecision(DEFAULT_PRECISION);

	gtStr
	<< ":";
	const HMM_PP::ModelParams::DataVal* meanRD = gt->getMeanRD();
	if (meanRD != NULL)
		gtStr
		<< *meanRD;
	else
		gtStr
		<< MISSING_VAL;

	if (hasOrigData()) {
		gtStr
		<< ":";
		const HMM_PP::ModelParams::DataVal* meanOrigRD = gt->getMeanOrigRD();
		if (meanOrigRD != NULL)
			gtStr
			<< *meanOrigRD;
		else
			gtStr
			<< MISSING_VAL;
	}

	char discoveredInSample = 'N';
	const set<string>* discoveredSamps = cnv->getDiscoveredSamples();
	if (discoveredSamps != NULL && discoveredSamps->find(gt->getSample()) != discoveredSamps->end())
		discoveredInSample = 'Y';
	gtStr
	<< ":" << discoveredInSample;

	return gtStr.str();
}

void XHMM::GenotypeOutputManager::initOstream(ostream& stream) {
	stream
	<< setiosflags(ios::fixed)
	<< setprecision(DEFAULT_PRECISION);
}
