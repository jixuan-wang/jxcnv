#include "GATKdepthMergerOutputManager.hpp"
#include "xhmmOutputManager.hpp"
#include "AuxiliaryCNVoutputManager.hpp"
#include "xhmmInputManager.hpp"
#include "OnDiskMatrixTransposer.hpp"

#include <Exception.hpp>
#include <NamedMatrix.hpp>
#include <utils.hpp>

#include <sstream>
#include <list>
#include <string>
using namespace std;

#define USE_DB_NO_STDIN "TMP_MAT_DB"

XHMM::GATKdepthMergerOutputManager::GATKdepthMergerOutputManager(list<string>* GATKdepthFiles, list<string>* GATKdepthListFiles, xhmmInputManager::OneToOneStrMap* mapSamples, string columnSuffix)
: _columnSuffix(columnSuffix), _allTargets(NULL), _samples(NULL), _readers(new list<GATKdepthReader*>()) {

	for (list<string>::const_iterator listIt = GATKdepthListFiles->begin(); listIt != GATKdepthListFiles->end(); ++listIt) {
		const string& listFile = *listIt;
		cerr << "Reading list of GATK depth files from " << listFile << endl;

		list<string>* filesToAdd = xhmmInputManager::readStringsFromFile(listFile);
		GATKdepthFiles->insert(GATKdepthFiles->end(), filesToAdd->begin(), filesToAdd->end());
		delete filesToAdd;
	}
	delete GATKdepthListFiles;

	if (GATKdepthFiles->empty())
		throw new Exception("Must provide at least one GATK depth file to merge");

	// Read the 1st GATK DoC file to get the "official" target list:
	GATKdepthReader* firstReader = new GATKdepthReader(*(GATKdepthFiles->begin()), columnSuffix, false);
	_allTargets = firstReader->getRemainingTargets();
	delete firstReader;

	list<string>* sampsList = new list<string>();

	for (list<string>::const_iterator fileIt = GATKdepthFiles->begin(); fileIt != GATKdepthFiles->end(); ++fileIt) {
		GATKdepthReader* reader = new GATKdepthReader(*fileIt, columnSuffix);
		_readers->push_back(reader);

		const list<GATKdepthReader::SampleIndex>* samps = reader->getSamples();
		for (list<GATKdepthReader::SampleIndex>::const_iterator it = samps->begin(); it != samps->end(); ++it) {
			string sample = it->first;
			if (mapSamples != NULL) {
				xhmmInputManager::OneToOneStrMap::iter findSampIt = mapSamples->find(sample);
				if (findSampIt != mapSamples->end()) {
					string newSample = findSampIt->second;
					if (newSample != sample) {
						cerr << "Mapping sample ID " << sample << " to ID " << newSample << endl;
						sample = newSample;
					}
				}
			}
			sampsList->push_back(sample);
		}
	}

	_samples = new vector<string>(sampsList->begin(), sampsList->end());
	delete sampsList;

	delete GATKdepthFiles;
	if (mapSamples != NULL)
		delete mapSamples;
}

XHMM::GATKdepthMergerOutputManager::~GATKdepthMergerOutputManager() {
	delete _allTargets;
	delete _samples;

	for (list<GATKdepthReader*>::const_iterator it = _readers->begin(); it != _readers->end(); ++it)
		delete *it;
	delete _readers;
}

void XHMM::GATKdepthMergerOutputManager::mergeGATKdepths(string outRdFile, bool outputSamplesByTargets, int rdPrecision, bool transposeInMemory) {
	cerr << "Writing GATK read-depth matrix of ";
	if (outputSamplesByTargets)
		cerr << _samples->size() << " samples by " << _allTargets->size() << " targets";
	else
		cerr << _allTargets->size() << " targets by " << _samples->size();
	cerr << " to " << outRdFile << endl;

	string matName = "GATK." + _columnSuffix;
	HMM_PP::BaseRealMat* rdMat = NULL;
	HMM_PP::ostreamWriter* outStream = NULL;
	OnDiskMatrixTransposer* matTransposer = NULL;

	if (!outputSamplesByTargets || !transposeInMemory) {
		outStream = HMM_PP::utils::getOstreamWriterFromFile(outRdFile);
		if (outStream == NULL)
			return;
		(*outStream)() << matName;
	}

	if (outputSamplesByTargets) {
		if (transposeInMemory) {
			rdMat = new HMM_PP::BaseRealMat(matName);
			rdMat->setDims(_samples->size(), _allTargets->size());

			for (uint sampInd = 0; sampInd < _samples->size(); ++sampInd)
				rdMat->setRowName(sampInd, (*_samples)[sampInd]);
		}
		else {
			// Will read in matrix (and cache to disk):
			matTransposer = new OnDiskMatrixTransposer((outRdFile != HMM_PP::utils::STD_STREAM ? outRdFile : USE_DB_NO_STDIN), _samples->size());
		}
	}
	else { // Targets x Samples:
		for (uint sampInd = 0; sampInd < _samples->size(); ++sampInd)
			(*outStream)() << '\t' << (*_samples)[sampInd];
		(*outStream)() << endl;

		(*outStream)() << setiosflags(ios::fixed);
		(*outStream)() << setprecision(rdPrecision);
	}

	uint targInd = 0;
	for (list<Interval>::const_iterator targIt = _allTargets->begin(); targIt != _allTargets->end(); ++targIt) {
		const Interval& targ = *targIt;
		if (outputSamplesByTargets) {
			if (transposeInMemory)
				rdMat->setColName(targInd, targ.intervalString());
			else
				(*outStream)() << '\t' << targIt->intervalString();
		}
		else // Targets x Samples:
			(*outStream)() << targIt->intervalString();

		uint sampInd = 0;
		list<string>* storeAllSampleVals = NULL;
		if (outputSamplesByTargets && !transposeInMemory)
			storeAllSampleVals = new list<string>();

		for (list<GATKdepthReader*>::const_iterator readIt = _readers->begin(); readIt != _readers->end(); ++readIt) {
			GATKdepthReader* reader = *readIt;
			if (!reader->hasNextTarget())
				throw new Exception("All GATK DoC files must have the same number of targets");

			GATKdepthReader::TargetSampleValues targSampVals = reader->nextTargetValues();
			if (targSampVals.first != targ)
				throw new Exception("All GATK DoC files must have the same targets");

			GATKdepthReader::SampleValues* sampVals = targSampVals.second;
			for (GATKdepthReader::SampleValues::const_iterator sampIt = sampVals->begin(); sampIt != sampVals->end(); ++sampIt) {
				if (outputSamplesByTargets) {
					if (transposeInMemory)
						(*rdMat)(sampInd++, targInd) = sampIt->second;
					else {
						stringstream valStr;
						valStr << setiosflags(ios::fixed);
						valStr << setprecision(rdPrecision);
						valStr << sampIt->second;

						storeAllSampleVals->push_back(valStr.str());
					}
				}
				else // Targets x Samples:
					(*outStream)() << '\t' << sampIt->second;
			}
			delete sampVals;
		}
		if (outputSamplesByTargets && !transposeInMemory)
			matTransposer->cacheNextRowToDBandDelete(storeAllSampleVals);
		if (!outputSamplesByTargets)
			(*outStream)() << endl;

		++targInd;
	}

	if (outputSamplesByTargets) {
		if (transposeInMemory)
			xhmmOutputManager::printMatrix(*rdMat, outRdFile, true, true, rdPrecision);
		else {
			(*outStream)() << endl;

			for (uint sampInd = 0; sampInd < _samples->size(); ++sampInd) {
				(*outStream)() << (*_samples)[sampInd];
				PrintOutputData* printSampData = new PrintOutputData((*outStream)());
				matTransposer->nextTransposedRowFromDB(printSampData);
				(*outStream)() << endl;
			}
		}
	}

	if (rdMat != NULL)
		delete rdMat;

	if (outStream != NULL)
		delete outStream;

	if (matTransposer != NULL)
		delete matTransposer;
}

/*
 * GATKdepthReader
 */
XHMM::GATKdepthMergerOutputManager::GATKdepthReader::GATKdepthReader(string GATKdepthFile, string columnSuffix, bool verbose)
: _stream(HMM_PP::utils::getIstreamLineReaderFromFile(GATKdepthFile)), _samples(new list<SampleIndex>()), _columnSuffix(columnSuffix), _columnSuffixLength(columnSuffix.size()) {

	if (verbose)
		cerr << "Reading GATK depth file " << GATKdepthFile << endl;

	readSamplesFromHeader();
	if (_samples->empty())
		throw new Exception("Could not read any samples from " + GATKdepthFile + " [expected suffix in header is: " + _columnSuffix + "]");
}

XHMM::GATKdepthMergerOutputManager::GATKdepthReader::~GATKdepthReader() {
	delete _stream;
	delete _samples;
}

void XHMM::GATKdepthMergerOutputManager::GATKdepthReader::readSamplesFromHeader() {
	if (_stream->eof())
		return;

	string* header = new string();
	*_stream >> *header;
	stringstream* headerStream = new stringstream(*header);
	delete header;

	uint ind = 0;
	while (*headerStream && !headerStream->eof()) {
		string field;
		// Instead of splitting by whitespace, use ONLY tab as delimiter for the header line (to allow for sample names with space in them):
		if (!getline(*headerStream, field, '\t'))
			break;

		if (ind != 0) { // since the first column contains the targets
			if (HMM_PP::endsWith(field, _columnSuffix)) {
				string sample = field.substr(0, field.size() - _columnSuffixLength);
				_samples->push_back(SampleIndex(sample, ind));
			}
			else if (HMM_PP::beginsWith(field, _columnSuffix)) { // Allow the "suffix" to be a prefix:
				string sample = field.substr(_columnSuffixLength);
				_samples->push_back(SampleIndex(sample, ind));
			}
		}

		++ind;
	}

	delete headerStream;
}

bool XHMM::GATKdepthMergerOutputManager::GATKdepthReader::hasNextTarget() {
	return !_stream->eof();
}

XHMM::GATKdepthMergerOutputManager::GATKdepthReader::TargetSampleValues XHMM::GATKdepthMergerOutputManager::GATKdepthReader::nextTargetValues(bool extractValues) {
	if (!hasNextTarget())
		throw new Exception("GATKdepthReader: Cannot get next if don't have any more");

	string* line = new string();
	*_stream >> *line;
	stringstream* lineStream = new stringstream(*line);
	delete line;

	uint ind = 0;
	string targetStr;
	*lineStream >> targetStr;
	if (!*lineStream)
		throw new Exception("GATKdepthReader: Could not read target");
	Interval targ(targetStr);

	SampleValues* sampVals = NULL;
	if (extractValues) {
		sampVals = new SampleValues();
		if (!_samples->empty()) {
			list<SampleIndex>::const_iterator sampIt = _samples->begin();

			while (*lineStream && !lineStream->eof()) {
				if (sampIt == _samples->end())
					break;

				if (++ind == sampIt->second) {
					string valString;
					*lineStream >> valString;
					if (!*lineStream)
						throw new Exception("GATKdepthReader: Could not read sample value");

					string useValString = valString;
					if (valString[0] == '>') // For example, if read-depth is ">5000", then parse out numeric value of "5000"
						useValString = valString.substr(1);
					stringstream valStr(useValString);
					double val;
					valStr >> val;
					if (!valStr)
						throw new Exception("GATKdepthReader: Could not read sample value: '" + valString + "'");
					sampVals->push_back(SampleValue(sampIt->first, val));

					++sampIt;
				}
				else {
					string dummy;
					*lineStream >> dummy;
					if (!*lineStream)
						throw new Exception("GATKdepthReader: Could not read next field");
				}
			}
		}
	}

	delete lineStream;

	return TargetSampleValues(targ, sampVals);
}

list<XHMM::Interval>* XHMM::GATKdepthMergerOutputManager::GATKdepthReader::getRemainingTargets() {
	list<Interval>* intervals = new list<Interval>();

	while (hasNextTarget())
		intervals->push_back(nextTargetValues(false).first);

	return intervals;
}
