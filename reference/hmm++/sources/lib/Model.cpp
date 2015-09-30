#include "Model.hpp"
#include "Simplex.hpp"
#include "Data.hpp"
#include "Exception.hpp"
#include "NamedVector.hpp"
#include "NamedMatrix.hpp"
#include "ModelParams.hpp"
#include "DataLoader.hpp"
#include "InferredDataFit.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

HMM_PP::Model::Model(ModelParams* modelParams)
: _modelParams(modelParams) {
}

HMM_PP::Model::~Model() {
	delete _modelParams;
}

void HMM_PP::Model::loadAllData() {
	while (_modelParams->hasNextLoadData()) {
		Data* d = _modelParams->loadNextData();
		cerr << "read " << d->getNumObservations() << " data-points for observation " << d->getId() << endl;
		addData(d);
	}
}

void HMM_PP::Model::_internal_setDataFitResult(InferredDataFit* idf, const uint ind) {
	this->setDataFitResult(idf, ind);
	finishDataAndFit(idf);
}

real HMM_PP::Model::fit(bool calcDigammas, bool viterbi, bool timeEvents) {
	real lik = 1;

	for (uint cnt = 0; cnt < getNumIndividuals(); cnt++) {
		const Data* seq = getData(cnt);
		pair<InferredDataFit*, real> inferLik = fitSequence(seq, calcDigammas, viterbi, timeEvents);
		lik *= inferLik.second;
		_internal_setDataFitResult(inferLik.first, cnt);
	}

	return lik;
}

pair<HMM_PP::InferredDataFit*, real> HMM_PP::Model::fitSequence(const Data* seq, bool calcDigammas, bool viterbi, bool timeEvents) {
	InferredDataFit* infer = new InferredDataFit(seq, this, timeEvents);

	cerr << "Running forward-backward for " << seq->getId() << endl;
	infer->forwardBackward(calcDigammas);

	if (viterbi) {
		cerr << "Running Viterbi for " << seq->getId() << endl;
		infer->viterbi();
	}

	return make_pair(infer, infer->calcLikelihood());
}

ostream& HMM_PP::Model::printModel(ostream& stream) const {
	stream << *_modelParams;
	return stream;
}

ostream& HMM_PP::Model::printExplicitModel(ostream& stream, const uint maxT) const {
	uint k = this->getNumHiddenStates();
	stream << "Initial probs:" << endl;
	for (uint i = 0; i < k; i++)
		stream << _modelParams->state(i) << '\t' << _modelParams->getInitProbs()[i] << endl;
	stream << endl;

	for (uint t = 0; t < maxT; ++t) {
		stream << "Transition matrix from t=" << t << " to t=" << t+1 << ":" << endl;
		NamedMatrix<real>* mat = transitionMatrix(t, t+1);
		stream << *mat << endl;
		delete mat;
	}
	return stream;
}

BaseReal HMM_PP::Model::transitionFunc(const uint i, const uint j, const uint t1, const uint t2) const {
	return _modelParams->transitionFunc(i, j, t1, t2);
}

HMM_PP::NamedMatrix<real>* HMM_PP::Model::transitionMatrix(const uint t1, const uint t2) const {
	NamedMatrix<real>* a = new NamedMatrix<real>();
	a->setDims(getNumHiddenStates(), getNumHiddenStates());

	for (uint i = 0; i < getNumHiddenStates(); ++i) {
		for (uint j = 0; j < getNumHiddenStates(); j++) {
			(*a)(i,j) = this->transitionFunc(i, j, t1, t2);
		}
	}

	return a;
}

HMM_PP::NamedMatrix<real>* HMM_PP::Model::transitionMatrixTranspose(const uint t1, const uint t2) const {
	NamedMatrix<real>* aTranspose = new NamedMatrix<real>();
	aTranspose->setDims(getNumHiddenStates(), getNumHiddenStates());

	for (uint i = 0; i < getNumHiddenStates(); ++i) {
		for (uint j = 0; j < getNumHiddenStates(); j++) {
			(*aTranspose)(i,j) = this->transitionFunc(j, i, t1, t2);
		}
	}

	return aTranspose;
}

real HMM_PP::Model::emissionFunc(const Data* seq, const uint hidden_state, const uint t) const {
	return _modelParams->emissionFunc(seq, hidden_state, t);
}

vector<real>* HMM_PP::Model::emissionFunc(const HMM_PP::Data* seq, const uint t) const {
	return _modelParams->emissionFunc(seq, t);
}

void HMM_PP::Model::baumWelch() {
	uint iter = 0;
	const uint mxiter = 500;
	const BaseReal EPS = 1e-8;
	real old_lik = 0;

	while (true)  {
		// E-step (considers all sequences)
		const bool CALC_DIGAMMAS = true; // need to calculate the digammas for updating the transition matrix
		const bool RUN_VITERBI = false;
		real lik = fit(CALC_DIGAMMAS, RUN_VITERBI);

		// M-step
		_modelParams->clear();
		reestimateTransitionAndEmissionNoScaling();
		_modelParams->scale(getNumIndividuals());

		// Re-estimate 'pi'
		reestimateInitProbs();

		// clear all data
		//clearAllData();

		// Next iteration?
		if (++iter == mxiter) break;

		cerr << "Baum-Welch iteration " << iter << '\t' << " LogLik = " << lik.getLog10Value() << endl;

		cerr << "PI";
		for (uint i = 0; i < _modelParams->getNumHiddenStates(); i++)
			cerr << '\t' << _modelParams->getInitProbs()[i];
		cerr << endl;

		cerr << *this;

		if (iter == 1 || (lik / old_lik).getLog10Value() > EPS)
			old_lik = lik;
		else
			break;

		cerr << "\n----------------------------------------------------------------\n\n";
	}
}

BaseReal f_simplex(HMM_PP::BaseRealVec& v, void* d) {
	HMM_PP::Model* model = (HMM_PP::Model*)d;

	model->getModelParams()->setVariableParams(v);
	real lik = model->fit(false, false); // don't need digammas or Viterbi
	model->reestimateInitProbs();

	// return the negative log-likelihood, where smaller numbers are "better" (higher likelihood of data)
	return - lik.getLog10Value();
}

void HMM_PP::Model::reestimateTransitionAndEmissionNoScaling() {
	// given an HMM, re-estimate parameters:
	for (uint cnt = 0; cnt < getNumIndividuals(); ++cnt) {
		InferredDataFit* fit = getDataFitResult(cnt);

		_modelParams->updateTransition(fit);
		_modelParams->updateEmission(fit);

		finishDataAndFit(fit);
	}
}

void HMM_PP::Model::reestimateInitProbs() {
	BaseRealVec initProbs(_modelParams->getNumHiddenStates(), 0);

	for (uint cnt = 0; cnt < getNumIndividuals(); ++cnt) {
		InferredDataFit* infer = getDataFitResult(cnt);

		for (uint i = 0; i< _modelParams->getNumHiddenStates(); i++)
			initProbs[i] += infer->getGamma(0, i).getValue();

		finishDataAndFit(infer);
	}

	_modelParams->setStartingProbsAndNormalize(initProbs);
}

bool HMM_PP::Model::nelderMead()  {
	Simplex simplex;

	simplex.set_function(f_simplex);
	simplex.set_data((void*)this);
	simplex.set_start(_modelParams->getVariableParams());
	simplex.optimize(1);

	return true;
}

void HMM_PP::Model::displayFittedSequence(ostream& stream) {
	for (uint cnt = 0; cnt < getNumIndividuals(); cnt++) {
		InferredDataFit* infer = getDataFitResult(cnt);

		infer->displayViterbiFittedSequence(stream);

		finishDataAndFit(infer);
	}
}
