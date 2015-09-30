#include "DiscoverOutputManager.hpp"
#include "SampleGenotypeQualsCalculator.hpp"
#include "CNVmodelParams.hpp"
#include "Genotype.hpp"
#include "AlleleQuals.hpp"
#include "AuxiliaryCNVoutputManager.hpp"
#include "ReadDepthMatrixLoader.hpp"
#include "PosteriorOutputManager.hpp"

#include <params.hpp>
#include <Model.hpp>
#include <Data.hpp>
#include <InferredDataFit.hpp>
#include <PreciseNonNegativeReal.hpp>

#include <set>
#include <fstream>
#include <iomanip>
#include <utility>
using namespace std;

const string XHMM::DiscoverOutputManager::SAMPLE = "SAMPLE";
const string XHMM::DiscoverOutputManager::CNV = "CNV";
const string XHMM::DiscoverOutputManager::INTERVAL = "INTERVAL";
const string XHMM::DiscoverOutputManager::EXACT_QUAL = "Q_EXACT";

XHMM::DiscoverOutputManager::DiscoverOutputManager
(string outFile, AuxiliaryCNVoutputManager* auxOutput,
		string posteriorFile,
		BaseReal discoverSomeQualThreshold, BaseReal maxScoreVal, int scoreDecimalPrecision,
		ReadDepthMatrixLoader* origDataLoader)
:   CNVoutputManager(origDataLoader),
    _outStream(HMM_PP::utils::getOstreamWriterFromFile(outFile)), _auxOutput(auxOutput),
    _posteriorFile(posteriorFile), _posteriorOutput(NULL),
    _discoverSomeQualThreshold(scoreDecimalPrecision, discoverSomeQualThreshold), _maxScoreVal(maxScoreVal), _scoreDecimalPrecision(scoreDecimalPrecision),
    _excludeSegmentTypes(new set<uint>()) {

	if (_maxScoreVal < discoverSomeQualThreshold)
		_maxScoreVal = discoverSomeQualThreshold;

	if (_outStream == NULL)
		throw new Exception("Cannot give empty file name for output");

	(*_outStream)() << setiosflags(ios::fixed);
	(*_outStream)() << setprecision(DEFAULT_PRECISION);

	printCNVheader();

	_excludeSegmentTypes->insert(CNVmodelParams::DIPLOID);
}

XHMM::DiscoverOutputManager::~DiscoverOutputManager() {
	delete _outStream;
	if (_posteriorOutput != NULL)
		delete _posteriorOutput;
	delete _excludeSegmentTypes;
}

void XHMM::DiscoverOutputManager::setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader) {
	CNVoutputManager::setModelAndReadDepthLoader(model, dataLoader);

	if (_posteriorFile != "")
		_posteriorOutput = new PosteriorOutputManager(_posteriorFile, CNVmodelParams::ALL_CN_TYPES, _model, _dataLoader);
}

void XHMM::DiscoverOutputManager::printAllCNVs(HMM_PP::Model* model) {
	for (uint i = 0; i < model->getNumIndividuals(); i++) {
		HMM_PP::InferredDataFit* idf = model->getDataFitResult(i);
		if (idf == NULL)
			continue;

		printDataType(idf);
		model->finishDataAndFit(idf);
	}
}

void XHMM::DiscoverOutputManager::printDataType(const HMM_PP::InferredDataFit* idf) {
	HMM_PP::Data* origData = CNVoutputManager::getOrigDataForNextDataFitOutput(idf);

	vector<HMM_PP::ullintPair> seqs = idf->getViterbiSegments(_excludeSegmentTypes, &(_dataLoader->getChrStopTargetIndices()));
	SampleGenotypeQualsCalculator* genotyper = new SampleGenotypeQualsCalculator(idf, _model, _maxScoreVal, _dataLoader);

	for (uint s = 0; s < seqs.size(); s++) {
		const uint t1 = seqs[s].first;
		const uint t2 = seqs[s].second;
		const uint type = idf->getViterbiPath(t1);

		printCNV(genotyper, t1, t2, type, origData);
	}

	if (origData != NULL)
		delete origData;

	delete genotyper;

	if (_posteriorOutput != NULL)
		_posteriorOutput->printSamplePosteriors(idf);
}

void XHMM::DiscoverOutputManager::printCNVheader() {
	ostream& outStream = (*_outStream)();

	outStream
	<< SAMPLE
	<< '\t' << CNV
	<< '\t' << INTERVAL
	<< '\t' << "KB"
	<< '\t' << "CHR"
	<< '\t' << "MID_BP"
	<< '\t' << "TARGETS"
	<< '\t' << "NUM_TARG"
	<< '\t' << EXACT_QUAL
	<< '\t' << "Q_SOME"
	<< '\t' << "Q_NON_DIPLOID"
	<< '\t' << "Q_START"
	<< '\t' << "Q_STOP"
	<< '\t' << "MEAN_RD";

	if (hasOrigData())
		outStream
		<< '\t' << "MEAN_ORIG_RD";

	outStream
	<< endl;
}

void XHMM::DiscoverOutputManager::printCNV(SampleGenotypeQualsCalculator* genotyper, const uint t1, const uint t2, const uint type, HMM_PP::Data* origData) {
	TargetRange range(t1, t2);
	Genotype* gt = new Genotype(genotyper, &range, CNVmodelParams::DIPLOID, type);
	const AlleleQuals* typeAllele = gt->getNonRefAllele(type);

	if (_discoverSomeQualThreshold.passScore(typeAllele->getHaveSomeEventScore())) {
		ostream& outStream = (*_outStream)();

		const Interval& mergedTarg = *(gt->getMergedTargets());
		const TargetRange& targRange = *(gt->getTargetRange());

		outStream
		<< gt->getSample()
		<< '\t' << _model->getModelParams()->state(typeAllele->getCNVtype())
		<< '\t' << mergedTarg
		<< '\t' << mergedTarg.spanKB()
		<< '\t' << mergedTarg.getChr()
		<< '\t' << mergedTarg.midpoint()
		<< '\t' << targRange.getOneBasedTargetIndexSpan()
		<< '\t' << targRange.getNumTargets()
		<< setprecision(_scoreDecimalPrecision)
		<< '\t' << typeAllele->getHaveExactEventScore()
		<< '\t' << typeAllele->getHaveSomeEventScore()
		<< '\t' << *(gt->getNonDiploidScore())
		<< '\t' << typeAllele->getStartStopScores().first
		<< '\t' << typeAllele->getStartStopScores().second
		<< setprecision(DEFAULT_PRECISION)
		<< '\t' << *(gt->getMeanRD());

		delete gt;

		if (origData != NULL) {
			const HMM_PP::ModelParams::DataVal* meanOrigRD = _model->calcRepresentativeDataVal(origData, t1, t2);
			outStream
			<< '\t' << *meanOrigRD;
			delete meanOrigRD;
		}

		outStream
		<< endl;

		if (_auxOutput != NULL)
			_auxOutput->printCNVtargets(genotyper, t1, t2, type, origData);
	}
}

XHMM::DiscoverOutputManager::IntervalToThreshMap* XHMM::DiscoverOutputManager::CNVtableToIntervals(xhmmInputManager::Table* cnvTable) {
	IntervalToThreshMap* targets = new IntervalToThreshMap();

	if (!cnvTable->hasColumn(INTERVAL)) {
		if (cnvTable->getNumColumns() != 1)
			throw new Exception("Input CNV table must have a column title " + INTERVAL);

		// Create a table with an 'INTERVAL' header if only had a single column:
		xhmmInputManager::Table::Row* colNames = new xhmmInputManager::Table::Row(1, INTERVAL);
		xhmmInputManager::Table* tmpCnvTable = new xhmmInputManager::Table(colNames);

		string firstInterval = cnvTable->getColumns()[0];
		tmpCnvTable->addRow(new xhmmInputManager::Table::Row(1, firstInterval));

		for (uint i = 0; i < cnvTable->getNumRows(); ++i) {
			string targString = cnvTable->getEntry(i, 0);
			tmpCnvTable->addRow(new xhmmInputManager::Table::Row(1, targString));
		}
		delete cnvTable;
		cnvTable = tmpCnvTable;
	}

	bool hasExactQuals = cnvTable->hasColumn(EXACT_QUAL);
	bool hasSampleInfo = cnvTable->hasColumn(SAMPLE);

	for (uint i = 0; i < cnvTable->getNumRows(); ++i) {
		string targString = cnvTable->getEntry(i, INTERVAL);
		Interval targ(targString);

		IntervalToThreshMap::iterator findIt = targets->find(targ);
		if (findIt == targets->end())
			findIt = targets->insert(make_pair(targ, ThreshAndSamples(NULL, NULL))).first;
		ThreshAndSamples& threshSamples = findIt->second;

		if (hasExactQuals) {
			if (threshSamples.first == NULL)
				threshSamples.first = new BaseReal(HMM_PP::PreciseNonNegativeReal<BaseReal>::REAL_INFINITY);
			BaseReal* thresh = threshSamples.first;

			stringstream qualStream;
			qualStream << cnvTable->getEntry(i, EXACT_QUAL);
			BaseReal qual;
			qualStream >> qual;
			if (!qualStream)
				throw new Exception("Unable to read quality from " + EXACT_QUAL + " field: " + qualStream.str());

			if (qual < *thresh)
				*thresh = qual;
		}

		if (hasSampleInfo) {
			if (threshSamples.second == NULL)
				threshSamples.second = new set<string>();
			threshSamples.second->insert( cnvTable->getEntry(i, SAMPLE) );
		}
	}

	delete cnvTable;

	return targets;
}
