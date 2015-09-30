#include "InferredDataFit.hpp"
#include "Model.hpp"
#include "ModelParams.hpp"
#include "Data.hpp"
#include "Timer.hpp"

#include <cmath>
#include <sstream>
#include <utility>
#include <algorithm>
using namespace std;

HMM_PP::InferredDataFit::InferredDataFit(const Data* data, Model* model, bool timeEvents)
: _data(data), _model(model),
  _alpha(NULL), _beta(NULL), _gamma(NULL), _digamma(NULL), _vpath(NULL),
  _runTimes(timeEvents ? new list<TimedEvent> : NULL) {
}

HMM_PP::InferredDataFit::~InferredDataFit() {
	if (_runTimes != NULL) {
		uint ind = 0;
		cerr << "Run-times for " << _data->getId() << ": ";
		for (list<TimedEvent>::const_iterator it = _runTimes->begin(); it != _runTimes->end(); ++it) {
			if (ind++ > 0)
				cerr << ", ";
			cerr << it->first << "= " << it->second;
		}
		cerr << endl;

		delete _runTimes;
	}

	delete _alpha;
	delete _beta;
	delete _gamma;
	delete _digamma;
	delete _vpath;
}

void HMM_PP::InferredDataFit::clear() {
	_alpha->clear();
	_beta->clear();
	_gamma->clear();
	_digamma->clear();

	_vpath->clear();
	vector<uint>().swap(*_vpath);

	if (_runTimes != NULL)
		_runTimes->clear();
}

/*
 * High-level HMM queries
 */
real HMM_PP::InferredDataFit::calcLikelihood() const {
	return calcLikelihood(*_alpha, *_beta, _data->getNumObservations() - 1);
}

real HMM_PP::InferredDataFit::calcLikelihood(const NamedMatrix<real>& alpha, const NamedMatrix<real>& beta, const uint t) const {
	const uint k = _model->getNumHiddenStates();

	real lk = 0;
	for (uint i = 0; i < k; i++)
		lk += alpha(t,i) * beta(t,i);

	return lk;
}

real HMM_PP::InferredDataFit::calcPosteriorNumerator(const uint t1, const vector<uint>& statePath) const {
	const uint t2 = _data->checkStartStopRange(t1, statePath);

	real posteriorNumerator = 1;

	posteriorNumerator *= (*_alpha)(t1, statePath[0]);
	posteriorNumerator *= (*_beta)(t2, statePath[statePath.size() - 1]);

	for (uint t = t1+1; t <= t2; ++t) {
		const uint lookupInd = t-t1;
		posteriorNumerator *= _model->transitionFunc(statePath[lookupInd-1], statePath[lookupInd], t-1, t);
		posteriorNumerator *= _model->emissionFunc(_data, statePath[lookupInd], t);
	}

	return posteriorNumerator;
}

real HMM_PP::InferredDataFit::calcPosterior(const uint t1, const vector<uint>& statePath) const {
	real posterior = calcPosteriorNumerator(t1, statePath);
	posterior /= calcLikelihood();
	return posterior;
}

real HMM_PP::InferredDataFit::calcPosteriorComplement(const uint t1, const vector<uint>& statePath) const {
	real posteriorNumerator = calcPosteriorNumerator(t1, statePath);
	real likelihood = calcLikelihood();

	real posteriorComplement = likelihood - posteriorNumerator;
	posteriorComplement /= likelihood;

	return posteriorComplement;
}

real HMM_PP::InferredDataFit::calcPosteriorNumerator(const uint t1, const uint t2, const uint state) const {
	_data->checkStartStopRange(t1, t2);

	vector<uint>* statePath = new vector<uint>(t2 - t1 + 1, state);
	real posteriorNumerator = calcPosteriorNumerator(t1, *statePath);
	delete statePath;

	return posteriorNumerator;
}

real HMM_PP::InferredDataFit::calcPosterior(const uint t1, const uint t2, const uint state) const {
	real posterior = calcPosteriorNumerator(t1, t2, state);
	posterior /= calcLikelihood();
	return posterior;
}

real HMM_PP::InferredDataFit::calcPosteriorComplement(const uint t1, const uint t2, const uint state) const {
	real posteriorNumerator = calcPosteriorNumerator(t1, t2, state);
	real likelihood = calcLikelihood();

	real posteriorComplement = likelihood - posteriorNumerator;
	posteriorComplement /= likelihood;

	return posteriorComplement;
}

/*
 * More sophisticated, high-level HMM queries
 */
real HMM_PP::InferredDataFit::calcLikelihoodGivenConstraints(const uint t1, const uint t2, Constraint* constraints) const {
	// Since doing forward-backward for HMM (i.e., "sum-product" message passing):
	NamedMatrix<real>::MarginalizeFunc margFunc = NamedMatrix<real>::sum;

	NamedMatrix<real>* alpha = new NamedMatrix<real>(t2 + 1, _model->getNumHiddenStates(), 0); // Ensure that (*alpha)[t2] can be used

	if (t1 > 0) // Ensure that the incoming forward message exists:
		(*alpha)[t1-1] = (*_alpha)[t1-1];

	for (uint t = t1; t <= t2; ++t)
		calcForwardMessage(t, *alpha, margFunc, constraints);
	delete constraints;

	real likelihoodGivenConstraints = calcLikelihood(*alpha, *_beta, t2);
	delete alpha;

	return likelihoodGivenConstraints;
}

// Calculates: P(y_1_T | x_t1 = statePath[0], ..., x_t2 = statePath[statePath.size()-1])
real HMM_PP::InferredDataFit::calcLikelihoodGivenEvent(const uint t1, const vector<uint>& statePath) const {
	const uint t2 = _data->checkStartStopRange(t1, statePath);

	TimeDependentConstraint* constraints = new TimeDependentConstraint();
	for (uint t = t1; t <= t2; ++t) {
		const uint lookupInd = t-t1;
		set<uint> permittedState;
		permittedState.insert(statePath[lookupInd]);
		constraints->addConstraintsAtTime(t, createConstraintsOnlyPermitStates(permittedState));
	}

	return calcLikelihoodGivenConstraints(t1, t2, constraints);
}

// Calculates: P(y_1_T | x_t1 = state, ..., x_t2 = state)
real HMM_PP::InferredDataFit::calcLikelihoodGivenEvent(const uint t1, const uint t2, const uint state) const {
	_data->checkStartStopRange(t1, t2);

	set<uint> permittedState;
	permittedState.insert(state);
	HomogeneousConstraint* constraints = createConstraintsOnlyPermitStates(permittedState);
	return calcLikelihoodGivenConstraints(t1, t2, constraints);
}

// Calculates numerator of: P(x_t \notin {excludeStates}, t1 <= t <= t2 | y_1_T)
real HMM_PP::InferredDataFit::calcPosteriorNumeratorAssignmentsExcludingStates(const uint t1, const uint t2, const set<uint>& excludeStates) const {
	_data->checkStartStopRange(t1, t2);

	HomogeneousConstraint* constraints = createConstraintsExcludedStates(excludeStates);
	return calcLikelihoodGivenConstraints(t1, t2, constraints);
}

real HMM_PP::InferredDataFit::calcPosteriorAssignmentsExcludingStates(const uint t1, const uint t2, const set<uint>& excludeStates) const {
	real posteriorAssignmentsExcludingStates = calcPosteriorNumeratorAssignmentsExcludingStates(t1, t2, excludeStates);
	posteriorAssignmentsExcludingStates /= calcLikelihood();
	return posteriorAssignmentsExcludingStates;
}

real HMM_PP::InferredDataFit::calcPosteriorComplementAssignmentsExcludingStates(const uint t1, const uint t2, const set<uint>& excludeStates) const {
	real posteriorNumeratorAssignmentsExcludingStates = calcPosteriorNumeratorAssignmentsExcludingStates(t1, t2, excludeStates);
	real likelihood = calcLikelihood();

	real posteriorComplement = likelihood - posteriorNumeratorAssignmentsExcludingStates;
	posteriorComplement /= likelihood;

	return posteriorComplement;
}

/*
 * Basic HMM operations
 */
void HMM_PP::InferredDataFit::forwardBackward(bool calcDigammas) {
	Timer* t = NULL;
	if (_runTimes != NULL)
		t = new Timer();

	if (_alpha != NULL)
		delete _alpha;
	_alpha = new NamedMatrix<real>();

	if (_beta != NULL)
		delete _beta;
	_beta = new NamedMatrix<real>();

	if (_gamma != NULL)
		delete _gamma;
	_gamma = new NamedMatrix<real>();

	forwardBackward_Viterbi(NamedMatrix<real>::sum, *_alpha, *_beta, *_gamma, true);

	if (calcDigammas)
		calcDigamma();

	if (_runTimes != NULL) {
		addCalcTime("FB", t->getDuration());
		delete t;
	}
}

void HMM_PP::InferredDataFit::viterbi() {
	Timer* t = NULL;
	if (_runTimes != NULL)
		t = new Timer();

	NamedMatrix<real>* alpha = new NamedMatrix<real>();
	NamedMatrix<real>* beta = new NamedMatrix<real>();
	NamedMatrix<real>* gamma = new NamedMatrix<real>();

	forwardBackward_Viterbi(NamedMatrix<real>::max, *alpha, *beta, *gamma, false);
	delete alpha;
	delete beta;

	// Retrieve the best path:
	if (_vpath != NULL)
		delete _vpath;
	_vpath = new vector<uint>(_data->getNumObservations());

	for (uint t = 0; t < _data->getNumObservations(); ++t)
		(*_vpath)[t] = (*gamma)[t].calcMaxArgmax().second;

	delete gamma;

	if (_runTimes != NULL) {
		addCalcTime("VT", t->getDuration());
		delete t;
	}
}

void HMM_PP::InferredDataFit::forwardBackward_Viterbi(NamedMatrix<real>::MarginalizeFunc margFunc, NamedMatrix<real>& alpha, NamedMatrix<real>& beta, NamedMatrix<real>& gamma, bool normalizeGamma) {
	const uint k = _model->getNumHiddenStates();
	alpha.setDims(_data->getNumObservations(), k, 0);
	beta.setDims(_data->getNumObservations(), k, 0);

	// compute forward probabilities at time t using the alpha values for time t-1:
	for (uint t = 0; t < _data->getNumObservations(); ++t) {
#if defined(DEBUG)
		cerr << "Calculating forward message " << t << endl;
#endif
		calcForwardMessage(t, alpha, margFunc);
#if defined(DEBUG)
		cerr << alpha[t] << endl;
#endif
	}

	// compute backward probabilities at time t using the beta values for time t+1:
	for (int t = _data->getNumObservations()-1; t >= 0; --t) {
#if defined(DEBUG)
		cerr << "Calculating backward message " << t << endl;
#endif
		calcBackwardMessage(t, beta, margFunc);
#if defined(DEBUG)
		cerr << beta[t] << endl;
#endif
	}

	calcGamma(alpha, beta, gamma, normalizeGamma);
}

void HMM_PP::InferredDataFit::calcGamma(const NamedMatrix<real>& alpha, const NamedMatrix<real>& beta, NamedMatrix<real>& gamma, bool normalize) {
	const uint k = _model->getNumHiddenStates();
	gamma.setDims(_data->getNumObservations(), k, 0);

	for (uint t = 0; t < _data->getNumObservations(); ++t) {
		gamma[t] = alpha[t];
		gamma[t] *= beta[t];
		if (normalize)
			gamma[t].normalize();
	}

#if defined(DEBUG)
	for (uint t = 0; t < _data->getNumObservations(); t++) {
		cerr << "a,b,g = " << t; // << '\t' << obs(t);
		for (uint i = 0; i < k; i++)
			cerr << '\t' << alpha(t,i) << " " << beta(t,i) << " " << gamma(t,i);
		cerr << endl;
	}
#endif
}

void HMM_PP::InferredDataFit::calcDigamma() {
	if (_digamma != NULL)
		delete _digamma;
	_digamma = new digamma_t(_data->getNumObservations());

	const uint k = _model->getNumHiddenStates();

	// calculate di-gammas given alphas and betas
	for (uint t = 0; t < _data->getNumObservations()-1; t++) {
		vector<real>* emissProbs = _model->emissionFunc(_data, t+1);
		real denom = 0;

		for (uint i = 0; i < k; i++)
			for (uint j = 0; j < k; j++)
				denom += (*_alpha)(t,i) * _model->transitionFunc(i,j,t,t+1) * (*emissProbs)[j] * (*_beta)(t+1,j);

		(*_digamma)[t].setDims(k, k);
		for (uint i = 0; i < k; i++) {
			for (uint j = 0; j < k; j++) {
				const real d = ((*_alpha)(t,i) * _model->transitionFunc(i,j,t,t+1) * (*emissProbs)[j] * (*_beta)(t+1,j)) / denom;
				(*_digamma)[t](i,j) = d;
			}
		}

		delete emissProbs;
	}
}

void HMM_PP::InferredDataFit::calcForwardMessage(const uint t, NamedMatrix<real>& alpha, NamedMatrix<real>::MarginalizeFunc margFunc, const Constraint* constraints) const {
	NamedVector<real>& forwardMess = alpha[t];

	if (t == 0) {
		for (uint i = 0; i < _model->getNumHiddenStates(); ++i)
			forwardMess[i] = _model->getModelParams()->getInitProbs()[i];
	}
	else {
		// TODO: an unnecessary copy of the incomingForwardMess vector here!
		const NamedMatrix<real>* incomingForwardMess = alpha[t-1].asMatrix();
		NamedMatrix<real>* transMat = _model->transitionMatrix(t-1, t);

		NamedMatrix<real>* matMultiply = NamedMatrix<real>::multiplyMarginalizeMatrices(*incomingForwardMess, *transMat, margFunc);
		delete incomingForwardMess;
		delete transMat;

		forwardMess = (*matMultiply)[0];
		delete matMultiply;
	}

	multiplyMessageTimesEmissionProbs(forwardMess, t, constraints);
}

void HMM_PP::InferredDataFit::calcBackwardMessage(const uint t, NamedMatrix<real>& beta, NamedMatrix<real>::MarginalizeFunc margFunc, const Constraint* constraints) const {
	NamedVector<real>& backwardMess = beta[t];

	if (t == _data->getNumObservations()-1) {
		backwardMess = 1.0;
	}
	else {
		NamedMatrix<real>* incomingBackwardMess = beta[t+1].asMatrix();
		multiplyMessageTimesEmissionProbs((*incomingBackwardMess)[0], t+1, constraints);

		NamedMatrix<real>* transMatTranspose = _model->transitionMatrixTranspose(t, t+1);

		NamedMatrix<real>* matMultiply = NamedMatrix<real>::multiplyMarginalizeMatrices(*incomingBackwardMess, *transMatTranspose, margFunc);
		delete incomingBackwardMess;
		delete transMatTranspose;

		backwardMess = (*matMultiply)[0];
		delete matMultiply;
	}
}

void HMM_PP::InferredDataFit::multiplyMessageTimesEmissionProbs(NamedVector<real>& message, const uint t, const Constraint* constraints) const {
	const uint k = _model->getNumHiddenStates();

	vector<real>* emissProbs = _model->emissionFunc(_data, t);

	for (uint i = 0; i < k; ++i) {
		real useEmissionProb;
		if (constraints == NULL || constraints->stateIsPermitted(i, t))
			useEmissionProb = (*emissProbs)[i];
		else
			useEmissionProb = 0;

		message[i] *= useEmissionProb;
	}

	delete emissProbs;
}

HMM_PP::InferredDataFit::HomogeneousConstraint* HMM_PP::InferredDataFit::createConstraintsExcludedStates(const set<uint>& excludedStates) const {
	return new HomogeneousConstraint(new set<uint>(excludedStates));
}

HMM_PP::InferredDataFit::HomogeneousConstraint* HMM_PP::InferredDataFit::createConstraintsOnlyPermitStates(const set<uint>& permittedStates) const {
	set<uint>* excludedStates = new set<uint>();
	insert_iterator<set<uint> > resultIt(*excludedStates, excludedStates->begin());

	set<uint>* allStates = getAllStatesList();
	set_difference(allStates->begin(), allStates->end(),
			permittedStates.begin(), permittedStates.end(),
			resultIt);
	delete allStates;

	return new HomogeneousConstraint(excludedStates);
}

set<uint>* HMM_PP::InferredDataFit::getAllStatesList() const {
	set<uint>* allStates = new set<uint>();

	const uint k = _model->getNumHiddenStates();
	for (uint i = 0; i < k; ++i)
		allStates->insert(i);

	return allStates;
}

real HMM_PP::InferredDataFit::getAlpha(const uint t, const uint i) const {
	return (*_alpha)(t,i);
}

real HMM_PP::InferredDataFit::getBeta(const uint t, const uint i) const {
	return (*_beta)(t,i);
}

real HMM_PP::InferredDataFit::getGamma(const uint t, const uint i) const {
	return (*_gamma)(t,i);
}

bool HMM_PP::InferredDataFit::hasDigamma() const {
	return _digamma != NULL && _digamma->size() == _data->getNumObservations();
}

real HMM_PP::InferredDataFit::getDigamma(const uint t, const uint i, const uint j) const {
	if (!hasDigamma())
		throw new Exception("Cannot get digamma results if not calculated");

	return (*_digamma)[t](i,j);
}

bool HMM_PP::InferredDataFit::hasViterbi() const {
	return _vpath != NULL && _vpath->size() == _data->getNumObservations();
}

uint HMM_PP::InferredDataFit::getViterbiPath(const uint t) const {
	if (!hasViterbi())
		throw new Exception("Cannot get Viterbi results if algorithm was not run");

	return (*_vpath)[t];
}

vector<HMM_PP::ullintPair> HMM_PP::InferredDataFit::getViterbiSegments(const set<uint>* excludeSegmentTypes, const set<uint>* forceSegmentBreakInds) const {
	if (!hasViterbi())
		throw new Exception("Cannot get Viterbi results if algorithm was not run");

	return calcSegments(*_vpath, excludeSegmentTypes, forceSegmentBreakInds);
}

void HMM_PP::InferredDataFit::insertSegments(vector<HMM_PP::ullintPair>& segments, const uint start, const uint stop, const set<uint>* forceSegmentBreakInds) {
	uint curStart = start;

	if (forceSegmentBreakInds != NULL) {
		for (set<uint>::const_iterator breakIt = forceSegmentBreakInds->lower_bound(start); breakIt != forceSegmentBreakInds->upper_bound(stop); ++breakIt) {
			const uint forceBreakInd = *breakIt;
			segments.push_back(ullintPair(curStart, forceBreakInd));

			curStart = forceBreakInd + 1;
		}
	}

	if (curStart <= stop)
		segments.push_back(ullintPair(curStart, stop));
}

void HMM_PP::InferredDataFit::displayViterbiFittedSequence(ostream& stream) const {
	const uint k = _model->getNumHiddenStates();

	for (uint t = 0; t < _data->getNumObservations(); t++) {
		stream << _data->getId() << '\t' << t << '\t';

		uint bestStateInd = 0;
		real best = -real::REAL_INFINITY;
		for (uint i = 1; i < k; i++) {
			if ((*_gamma)(t,i) > best) {
				best = (*_gamma)(t,i);
				bestStateInd = i;
			}
		}
		string mstlik = _model->getModelParams()->state(bestStateInd);

		_data->printDatapoint(stream, t);

		if (hasViterbi())
			stream << '\t' << _model->getModelParams()->state(getViterbiPath(t));

		for (uint i = 0; i < k; i++)
			stream << '\t' << (*_gamma)(t,i);

		stream << endl;
	}
}

void HMM_PP::InferredDataFit::addCalcTime(const string calc, const uint timeInSec) const {
	if (isCalcTimes())
		_runTimes->push_back(TimedEvent(calc, timeInSec));
}
