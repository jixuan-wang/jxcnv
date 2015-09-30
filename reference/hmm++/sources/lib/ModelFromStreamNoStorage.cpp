#include "ModelFromStreamNoStorage.hpp"
#include "Data.hpp"
#include "Exception.hpp"
#include "DataOutputter.hpp"
#include "InferredDataFit.hpp"

HMM_PP::ModelFromStreamNoStorage::ModelFromStreamNoStorage(ModelParams* modelParams, DataOutputter<InferredDataFit>* dataOuputter)
: Model(modelParams), _dataOuputter(dataOuputter) {
}

HMM_PP::ModelFromStreamNoStorage::~ModelFromStreamNoStorage() {
}

real HMM_PP::ModelFromStreamNoStorage::fit(bool calcDigammas, bool viterbi, bool timeEvents) {
	real lik = 1;

	while (_modelParams->hasNextLoadData()) {
		Data* seq = _modelParams->loadNextData();

		pair<InferredDataFit*, real> inferLik = fitSequence(seq, calcDigammas, viterbi, timeEvents);
		lik *= inferLik.second;

		InferredDataFit* fit = inferLik.first;
		_dataOuputter->printDataType(fit);

		delete fit;
		delete seq;
	}

	return lik;
}
