#include "CNVoutputManager.hpp"
#include "ReadDepthMatrixLoader.hpp"

#include <Model.hpp>
#include <Data.hpp>
#include <InferredDataFit.hpp>

#include <sstream>
using namespace std;

XHMM::CNVoutputManager::CNVoutputManager(ReadDepthMatrixLoader* origDataLoader)
: _model(NULL), _dataLoader(NULL),
  _origDataLoader(origDataLoader) {
}

XHMM::CNVoutputManager::~CNVoutputManager() {
}

void XHMM::CNVoutputManager::setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader) {
	if (_model != NULL || _dataLoader != NULL)
		throw new Exception("Cannot call setModelAndReadDepthLoader() more than once");
	_model = model;
	_dataLoader = dataLoader;

	if (hasOrigData()) {
		if (_origDataLoader->getNumTargets() != _dataLoader->getNumTargets()) {
			stringstream str;
			str << "Must provide original read depth data matrix of identical size - expected: " << _dataLoader->getNumTargets() << " found: " << _origDataLoader->getNumTargets();
			throw new Exception(str.str());
		}

		for (uint i = 0; i < _dataLoader->getNumTargets(); ++i) {
			if (_origDataLoader->getTarget(i) != _dataLoader->getTarget(i)) {
				stringstream str;
				str << "Must provide original read depth data matrix of identical targets - expected: " << _dataLoader->getTarget(i) << " found: " << _origDataLoader->getTarget(i);
				throw new Exception(str.str());
			}
		}
	}
}

HMM_PP::Data* XHMM::CNVoutputManager::getOrigDataForNextDataFitOutput(const HMM_PP::InferredDataFit* idf) {
	if (_model == NULL || _dataLoader == NULL)
		throw new Exception("Must call setModelAndReadDepthLoader() before using object");

	HMM_PP::Data* origData = NULL;
	if (hasOrigData()) {
		if (!_origDataLoader->hasNext())
			throw new Exception("Must provide original read depth data matrix of identical samples");
		origData = _origDataLoader->next();

		if (origData->getId() != idf->getData()->getId()) {
			stringstream str;
			str << "Must provide original read depth data matrix of samples in identical order. Expected: " << idf->getData()->getId() << ", got " << origData->getId();
			throw new Exception(str.str());
		}
	}

	return origData;
}
