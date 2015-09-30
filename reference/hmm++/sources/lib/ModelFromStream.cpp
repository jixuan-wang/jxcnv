#include "ModelFromStream.hpp"
#include "Data.hpp"
#include "DataLoader.hpp"
#include "Exception.hpp"
#include "InferredDataFit.hpp"

#include <sstream>
using namespace std;

HMM_PP::ModelFromStream::ModelFromStream(ModelParams* modelParams)
: Model(modelParams), _sequences() {
	this->loadAllData();
}

HMM_PP::ModelFromStream::~ModelFromStream() {
	for (uint s = 0; s < _sequences.size(); ++s) {
		const DataAndFit& df = _sequences[s];
		Data* seq = df.first;
		InferredDataFit* fit = df.second;

		if (seq != NULL)
			delete seq;
		if (fit != NULL)
			delete fit;
	}
	_sequences.clear();
}

void HMM_PP::ModelFromStream::loadAllData() {
	_sequences.clear();
	HMM_PP::Model::loadAllData();
}

uint HMM_PP::ModelFromStream::addData(HMM_PP::Data* d) {
	_sequences.push_back(DataAndFit(d, NULL));
	return _sequences.size() - 1;
}

const HMM_PP::Data* HMM_PP::ModelFromStream::getData(const uint i) {
	if (i >= _sequences.size()) {
		stringstream str;
		str << "Cannot request non-existent individual " << i;
		throw new Exception(str.str());
	}

	return _sequences[i].first;
}

HMM_PP::InferredDataFit* HMM_PP::ModelFromStream::getDataFitResult(const uint i) {
	if (i >= _sequences.size()) {
		stringstream str;
		str << "Cannot request non-existent individual " << i;
		throw new Exception(str.str());
	}

	InferredDataFit* fit = _sequences[i].second;
	if (fit == NULL) {
		stringstream str;
		str << "Cannot request non-existent fit results for individual " << i;
		throw new Exception(str.str());
	}

	return fit;
}

void HMM_PP::ModelFromStream::setDataFitResult(InferredDataFit* idf, const uint ind) {
	if (ind >= _sequences.size() || idf == NULL) {
		stringstream str;
		str << "Cannot set fit results for individual " << ind;
		throw new Exception(str.str());
	}

	if (_sequences[ind].second != NULL)
		delete _sequences[ind].second;

	_sequences[ind].second = idf;
}
