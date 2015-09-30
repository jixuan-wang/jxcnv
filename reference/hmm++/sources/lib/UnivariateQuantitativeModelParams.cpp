#include "UnivariateQuantitativeModelParams.hpp"
#include "utils.hpp"
#include "Model.hpp"
#include "Data.hpp"
#include "InferredDataFit.hpp"

HMM_PP::UnivariateQuantitativeModelParams::UnivariateQuantitativeModelParams(istream& stream, HMM_PP::DataLoader<UnivariateQuantitativeData>* dataLoader)
: ModelParams(new HiddenStateParams(stream)),
  HomogeneousModelParams(ModelParams::_params),
  HMM_PP::ModelParamsData<UnivariateQuantitativeData>(ModelParams::_params, dataLoader),
  _mean(), _var() {

	constructVariableParamsFromStream(stream);
	_params->setStartingProbsAndNormalize(this->getHomogeneousStationaryProbs());
}

void HMM_PP::UnivariateQuantitativeModelParams::constructVariableParamsFromStream(istream& stream) {
	BaseRealVec paramValues;

	HomogeneousModelParams::readTransMatParams(stream, _params->getNumHiddenStates(), paramValues);
	readGaussianParams(stream, _params->getNumHiddenStates(), paramValues);

	setVariableParams(paramValues);
}

void HMM_PP::UnivariateQuantitativeModelParams::readGaussianParams(istream& stream, const uint numHiddenStates, BaseRealVec& paramValues) {
	// get means and sd for each class as starting values/fixed parameters
	for (uint j = 0; j < numHiddenStates; j++) {
		BaseReal mean, sd;
		stream >> mean >> sd;
		if (!stream)
			throw new Exception("Failed to read mean and variance from stream");
		paramValues.addElement(mean);
		paramValues.addElement(sd);
	}
}

void HMM_PP::UnivariateQuantitativeModelParams::getGaussianParams(BaseRealVec& v) const {
	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		v.addElement(_mean[i]);
		v.addElement(sqrt(_var[i]));
	}
}

void HMM_PP::UnivariateQuantitativeModelParams::setGaussianFromParams(const BaseRealVec& v, uint& vInd) {
	// expecting under saturated model, k*2 for 'b':
	const uint remainingSize = (v.size() - vInd);
	if (remainingSize < 2 * _params->getNumHiddenStates()) {
		stringstream str;
		str << "Internal error in UnivariateQuantitativeModelParams::setVariableParams()" << remainingSize << " " << _params->getNumHiddenStates() << endl;
		throw new Exception(str.str());
	}
	_mean.clear();
	_mean.resize(_params->getNumHiddenStates(), 0);

	_var.clear();
	_var.resize(_params->getNumHiddenStates(), 1);

	// set mean/var for emission probs [means/SDs are not constrained]:
	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		_mean[i] = v[vInd++];
		BaseReal sd = v[vInd++];
		_var[i] = sd * sd;

		if (_var[i] <= 0)
			_var[i] = PreciseNonNegativeReal<BaseReal>::REAL_EPSILON;
	}
}

HMM_PP::BaseRealVec HMM_PP::UnivariateQuantitativeModelParams::getVariableParams() const {
	BaseRealVec v;

	getTransMatParams(v);
	getGaussianParams(v);

	return v;
}

void HMM_PP::UnivariateQuantitativeModelParams::setVariableParams(const BaseRealVec& v) {
	uint vInd = 0;
	setTransMatFromParams(v, vInd);
	setGaussianFromParams(v, vInd);
}

#define EMISS_FUNC dnorm(d->val(t), _mean[i], _var[i])

real HMM_PP::UnivariateQuantitativeModelParams::emissionFunc(const HMM_PP::Data* sequence, const uint i, const uint t) const {
	const DataType* d = castData(sequence);
	return EMISS_FUNC;
}

vector<real>* HMM_PP::UnivariateQuantitativeModelParams::emissionFunc(const HMM_PP::Data* sequence, const uint t) const {
	const DataType* d = castData(sequence);
	vector<real>* emissProbs = new vector<real>(getNumHiddenStates());

	for (uint i = 0; i < getNumHiddenStates(); ++i)
		(*emissProbs)[i] = EMISS_FUNC;

	return emissProbs;
}

void HMM_PP::UnivariateQuantitativeModelParams::clear() {
	HomogeneousModelParams::clear();
	_mean = 0;
	_var = 0;
}

void HMM_PP::UnivariateQuantitativeModelParams::updateEmission(const HMM_PP::InferredDataFit* infer) {
	// estimate mean and variance | class
	const DataType* d = castData(infer->getData());
	const uint n = d->getNumObservations();

	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		real numer = 0;
		real denom = 0;
		for (uint t = 0; t < n; t++) {
			numer += real(d->val(t)) * infer->getGamma(t,i);
			denom += infer->getGamma(t,i);
		}
		_mean[i] = (numer / denom).getValue();

		numer = denom = 0;
		for (uint t = 0; t < n; t++) {
			numer += infer->getGamma(t,i) * real((d->val(t) - _mean[i]) * (d->val(t) - _mean[i]));
			denom += infer->getGamma(t,i);
		}

		_var[i] = (numer / denom).getValue();
	}
}

void HMM_PP::UnivariateQuantitativeModelParams::scale(const BaseReal n) {
	HomogeneousModelParams::scale(n);
	_mean /= n;
	_var /= (n * n);
}

ostream& HMM_PP::UnivariateQuantitativeModelParams::print(ostream& stream) const {
	HMM_PP::HomogeneousModelParams::print(stream)
	<< "U" << endl
	<< _mean << endl
	<< "S" << endl
	<< _var << endl;
	return stream;
}

HMM_PP::UnivariateQuantitativeModelParams::QuantitativeDataVal* HMM_PP::UnivariateQuantitativeModelParams::calcRepresentativeDataVal(const HMM_PP::Data* d, const uint t1, const uint t2) const {
	d->checkStartStopRange(t1, t2);
	const DataType* dt = castData(d);

	BaseReal sumOfRD = 0;
	uint numTargets = 0;
	for (uint t = t1; t <= t2; ++t) {
		sumOfRD += dt->val(t);
		++numTargets;
	}

	return new QuantitativeDataVal(sumOfRD / numTargets);
}
