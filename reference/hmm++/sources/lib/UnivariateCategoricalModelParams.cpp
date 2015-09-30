#include "UnivariateCategoricalModelParams.hpp"
#include "utils.hpp"
#include "Model.hpp"
#include "Data.hpp"
#include "Simplex.hpp"
#include "Exception.hpp"
#include "InferredDataFit.hpp"

HMM_PP::UnivariateCategoricalModelParams::UnivariateCategoricalModelParams(istream& stream, HMM_PP::DataLoader<UnivariateCategoricalData>* dataLoader)
: ModelParams(new HiddenStateParams(stream)),
  HomogeneousModelParams(ModelParams::_params),
  ModelParamsData<UnivariateCategoricalData>(ModelParams::_params, dataLoader) {

	constructFixedParamsFromStream(stream);
	constructVariableParamsFromStream(stream);
	_params->setStartingProbsAndNormalize(this->getHomogeneousStationaryProbs());
}

void HMM_PP::UnivariateCategoricalModelParams::constructFixedParamsFromStream(istream& stream) {
	stream >> _numObservedStates;
	if (!stream)
		throw new Exception("Failed to read number of observed categories from input stream");

	HMM_PP::CategoricalData::clearCategories();
	for (uint j = 0; j < _numObservedStates; j++) {
		string label;
		stream >> label;
		if (!stream)
			throw new Exception("Failed to read observed category names from input stream");
		HMM_PP::CategoricalData::ensureCategoryToType(label);
	}
	if (HMM_PP::CategoricalData::getNumCategories() != _numObservedStates) {
		stringstream str;
		str << "Read " << _numObservedStates << " categories, but only " << HMM_PP::CategoricalData::getNumCategories() << " are unique";
		throw new Exception(str.str());
	}
}

void HMM_PP::UnivariateCategoricalModelParams::constructVariableParamsFromStream(istream& stream) {
	BaseRealVec paramValues;

	HomogeneousModelParams::readTransMatParams(stream, _params->getNumHiddenStates(), paramValues);

	// initial starting values for emission matrix [ordered by hidden state]:
	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		for (uint j = 0; j < _numObservedStates-1; j++) {
			BaseReal d;
			stream >> d;
			if (!stream)
				throw new Exception("Failed to read emission matrix from stream");
			paramValues.addElement(d);
		}

		// SKIP last entry of each row, since must normalize to a sum of 1 anyway:
		BaseReal dummy;
		stream >> dummy;
		if (!stream)
			throw new Exception("Failed to read emission matrix from stream");
	}

	setVariableParams(paramValues);
}

HMM_PP::BaseRealVec HMM_PP::UnivariateCategoricalModelParams::getVariableParams() const {
	BaseRealVec v;

	getTransMatParams(v);

	for (uint i = 0; i < _params->getNumHiddenStates(); i++)
		for (uint j = 0; j < _numObservedStates-1; j++)
			v.addElement(_emissionMat(i,j).getValue());

	return v;
}

void HMM_PP::UnivariateCategoricalModelParams::setVariableParams(const BaseRealVec& v) {
	uint vInd = 0;
	setTransMatFromParams(v, vInd);

	// expecting under saturated model, k*(s-1) for 'b':
	const uint remainingSize = (v.size() - vInd);
	if (remainingSize != _params->getNumHiddenStates() * (_numObservedStates - 1)) {
		stringstream str;
		str << "Internal error in UnivariateCategoricalModelParams::setVariableParams()" << remainingSize << " " << _params->getNumHiddenStates() << " " << _numObservedStates << endl;
		throw new Exception(str.str());
	}
	_emissionMat.setDims(_params->getNumHiddenStates(), _numObservedStates, 1/(BaseReal)_numObservedStates);

	// Read the emission matrix:
	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		BaseReal sum = 0;
		for (uint j = 0; j < _numObservedStates-1; j++) {
			_emissionMat(i,j) = v[vInd];
			sum += v[vInd];
			++vInd;
		}

		// constrain row to sum to 1.00
		if (sum > 1)  {
			for (uint j = 0; j < _numObservedStates-1; j++)
				_emissionMat(i,j) /= sum;
			sum = 1;
		}
		_emissionMat(i, _numObservedStates-1) = 1 - sum;
	}
}

#define EMISS_FUNC _emissionMat(i, d->val(t))

real HMM_PP::UnivariateCategoricalModelParams::emissionFunc(const HMM_PP::Data* sequence, const uint i, const uint t) const {
	const DataType* d = castData(sequence);
	return EMISS_FUNC;
}

vector<real>* HMM_PP::UnivariateCategoricalModelParams::emissionFunc(const HMM_PP::Data* sequence, const uint t) const {
	const DataType* d = castData(sequence);
	vector<real>* emissProbs = new vector<real>(getNumHiddenStates());

	for (uint i = 0; i < getNumHiddenStates(); ++i)
		(*emissProbs)[i] = EMISS_FUNC;

	return emissProbs;
}

void HMM_PP::UnivariateCategoricalModelParams::clear() {
	HomogeneousModelParams::clear();
	_emissionMat = 0;
}

void HMM_PP::UnivariateCategoricalModelParams::updateEmission(const HMM_PP::InferredDataFit* infer) {
	const DataType* d = castData(infer->getData());
	const uint n = d->getNumObservations();

	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		real denom = 0;
		for (uint t = 0; t < n; t++)
			denom += infer->getGamma(t,i);

		for (uint j = 0; j < _numObservedStates; j++) {
			real numer = 0;
			for (uint t = 0; t < n; t++) {
				if (d->val(t) == j)
					numer += infer->getGamma(t,i);
			}
			_emissionMat(i,j) += (numer / denom).getValue();
		}
	}
}

void HMM_PP::UnivariateCategoricalModelParams::scale(const BaseReal n) {
	HomogeneousModelParams::scale(n);
	_emissionMat /= n;
}

ostream& HMM_PP::UnivariateCategoricalModelParams::print(ostream& stream) const {
	HMM_PP::HomogeneousModelParams::print(stream)
	<< "B" << endl
	<< _emissionMat << endl;
	return stream;
}

HMM_PP::UnivariateCategoricalModelParams::CategoricalDataVal* HMM_PP::UnivariateCategoricalModelParams::calcRepresentativeDataVal(const HMM_PP::Data* d, const uint t1, const uint t2) const {
	d->checkStartStopRange(t1, t2);
	const DataType* dt = castData(d);

	map<uint, uint> catCounts;
	for (uint t = t1; t <= t2; ++t) {
		uint cat = dt->val(t);
		catCounts[cat]++;
	}

	uint argMaxCat = UINT_INFINITY;
	uint max = 0;
	for (map<uint, uint>::const_iterator it = catCounts.begin(); it != catCounts.end(); ++it) {
		if (it->second > max) {
			max = it->second;
			argMaxCat = it->first;
		}
	}

	return new CategoricalDataVal(argMaxCat);
}
