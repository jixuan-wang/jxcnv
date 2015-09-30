#include "CopyNumberVariant.hpp"
#include "AlleleQuals.hpp"
#include "ReadDepthMatrixLoader.hpp"

#include <sstream>
using namespace std;

XHMM::CopyNumberVariant::CopyNumberVariant(const vector<uint>& nonRefTypes, const HMM_PP::ullintPair* t1_t2, const ReadDepthMatrixLoader* dataLoader, Genotype::RealThresh genotypingThreshold, bool storeGenotypes, set<string>* discoveredSamples)
: _mergedTargets(t1_t2 == NULL ? NULL : new Interval(dataLoader->getTarget(t1_t2->first) + dataLoader->getTarget(t1_t2->second))),
  _prevTarget(NULL), _postTarget(NULL),
  _targRange(t1_t2 == NULL ? NULL : new TargetRange(*t1_t2)),
  _nonRefTypeToOrderIndex(new map<uint, uint>()),
  _sampleGenotypes(storeGenotypes ? new map<string, Genotype*>() : NULL),
  _nonRefTypeToSampleCount(new map<uint, uint>()), _calledSampleCount(0),
  _genotypingThreshold(genotypingThreshold),
  _discoveredSamples(discoveredSamples) {

	if (t1_t2 != NULL) {
		if (t1_t2->first > 0)
			_prevTarget = new Interval(dataLoader->getTarget(t1_t2->first - 1));
		if (t1_t2->second + 1 < dataLoader->getNumTargets())
			_postTarget = new Interval(dataLoader->getTarget(t1_t2->second + 1));
	}

	for (uint i = 0; i < nonRefTypes.size(); ++i)
		(*_nonRefTypeToOrderIndex)[nonRefTypes[i]] = i;
}

XHMM::CopyNumberVariant::~CopyNumberVariant() {
	if (_mergedTargets != NULL)
		delete _mergedTargets;

	if (_prevTarget != NULL)
		delete _prevTarget;
	if (_postTarget != NULL)
		delete _postTarget;

	if (_targRange != NULL)
		delete _targRange;

	delete _nonRefTypeToOrderIndex;

	if (_sampleGenotypes != NULL) {
		for (map<string, Genotype*>::const_iterator it = _sampleGenotypes->begin(); it != _sampleGenotypes->end(); ++it)
			delete it->second;
		delete _sampleGenotypes;
	}

	delete _nonRefTypeToSampleCount;

	if (_discoveredSamples != NULL)
		delete _discoveredSamples;
}

// NOTE: No error-checking as to the Genotype added:
void XHMM::CopyNumberVariant::addGenotype(string sample, Genotype* gt) {
	if (!gt->isNoCall())
		++_calledSampleCount;

	if (gt->isNonReference())
		++((*_nonRefTypeToSampleCount)[gt->getNonRefCall()->getCNVtype()]);

	if (_sampleGenotypes != NULL) {
		if (_sampleGenotypes->find(sample) != _sampleGenotypes->end())
			throw new Exception("Cannot add multiple genotypes for sample " + sample);
		(*_sampleGenotypes)[sample] = gt;
	}
	else
		delete gt;
}

const XHMM::Genotype* XHMM::CopyNumberVariant::getGenotype(string sample) const {
	if (_sampleGenotypes == NULL)
		throw new Exception("Cannot get genotype for sample: " + sample + " since not storing genotypes");

	map<string, Genotype*>::const_iterator findIt = _sampleGenotypes->find(sample);
	if (findIt == _sampleGenotypes->end())
		throw new Exception("Cannot get genotype for unknown sample: " + sample);

	return findIt->second;
}

uint XHMM::CopyNumberVariant::getNonRefTypeIndex(const uint type) const {
	map<uint, uint>::const_iterator findIt = _nonRefTypeToOrderIndex->find(type);
	if (findIt == _nonRefTypeToOrderIndex->end()) {
		stringstream str;
		str << "Cannot get index of invalid type: " << type;
		throw new Exception(str.str());
	}

	return findIt->second;
}

uint XHMM::CopyNumberVariant::getNonRefTypeCount(const uint type) const {
	map<uint, uint>::const_iterator findIt = _nonRefTypeToSampleCount->find(type);
	if (findIt == _nonRefTypeToSampleCount->end())
		return 0;

	return findIt->second;
}
