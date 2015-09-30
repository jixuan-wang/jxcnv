#include "HomogeneousModelParams.hpp"
#include "utils.hpp"
#include "Model.hpp"
#include "Data.hpp"
#include "Simplex.hpp"
#include "Exception.hpp"
#include "MatrixDecomp.hpp"
#include "PreciseNonNegativeReal.hpp"
#include "InferredDataFit.hpp"

HMM_PP::HomogeneousModelParams::HomogeneousModelParams(HiddenStateParams* params)
: ModelParams(params), _transitionMat() {
}

HMM_PP::HomogeneousModelParams::~HomogeneousModelParams() {
}

HMM_PP::BaseRealVec HMM_PP::HomogeneousModelParams::getHomogeneousStationaryProbs() const {
	HMM_PP::BaseRealMat* transMatTranspose = _transitionMat.transpose();
	HMM_PP::MatrixDecomp::EigenDecomp<BaseReal> decomp = HMM_PP::MatrixDecomp::rightEigen<BaseReal>(transMatTranspose);
	if (!PreciseNonNegativeReal<BaseReal>::equalsReal((*decomp.eigenValsReal)[0], 1.0) || !PreciseNonNegativeReal<BaseReal>::equalsReal((*decomp.eigenValsImaginary)[0], 0.0))
		throw new Exception("Transition matrix must be stochastic [hence with largest eigenvalue of 1]");

	HMM_PP::BaseRealMat* eigenVectorsAsRows = decomp.eigenVectors->transpose();
	decomp.deleteData();
	HMM_PP::BaseRealVec initProbs((*eigenVectorsAsRows)[0]);
	delete eigenVectorsAsRows;

	return initProbs.normalize();
}

void HMM_PP::HomogeneousModelParams::readTransMatParams(istream& stream, const uint numHiddenStates, BaseRealVec& paramValues) {
	for (uint i = 0; i < numHiddenStates; i++) {
		for (uint j = 0; j < numHiddenStates-1; j++) {
			BaseReal d;
			stream >> d;
			if (!stream)
				throw new Exception("Failed to read transition matrix from input stream");
			paramValues.addElement(d);
		}

		// SKIP last entry, since must normalize to a sum of 1 anyway:
		BaseReal dummy;
		stream >> dummy;
		if (!stream)
			throw new Exception("Failed to read transition matrix from input stream");
	}
}

void HMM_PP::HomogeneousModelParams::getTransMatParams(BaseRealVec& v) const {
	for (uint i = 0; i < _params->getNumHiddenStates(); i++)
		for (uint j = 0; j < _params->getNumHiddenStates()-1; j++)
			v.addElement(_transitionMat(i,j));
}

void HMM_PP::HomogeneousModelParams::setTransMatFromParams(const BaseRealVec& v, uint& vInd) {
	// expecting under saturated model, k*(k-1) for 'a':
	const uint remainingSize = (v.size() - vInd);
	if (remainingSize < _params->getNumHiddenStates() * (_params->getNumHiddenStates() - 1)) {
		stringstream str;
		str << "Internal error in HomogeneousModelParams::setTransMatFromParams()" << remainingSize << " " << _params->getNumHiddenStates() << endl;
		throw new Exception(str.str());
	}
	_transitionMat.setDims(_params->getNumHiddenStates(), _params->getNumHiddenStates(), 1/(BaseReal)_params->getNumHiddenStates());

	// Read the transition matrix:
	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		BaseReal sum = 0;
		for (uint j = 0; j < _params->getNumHiddenStates()-1; j++) {
			_transitionMat(i,j) = v[vInd];
			sum += v[vInd];
			++vInd;
		}

		// constrain row to sum to 1.00
		if (sum > 1)  {
			for (uint j = 0; j < _params->getNumHiddenStates()-1; j++)
				_transitionMat(i,j) /= sum;
			sum = 1;
		}
		_transitionMat(i, _params->getNumHiddenStates()-1) = 1 - sum;
	}
}

BaseReal HMM_PP::HomogeneousModelParams::transitionFunc(const uint i, const uint j, const uint t1, const uint t2) const {
	return _transitionMat(i,j);
}

void HMM_PP::HomogeneousModelParams::clear() {
	_transitionMat = 0;
}

void HMM_PP::HomogeneousModelParams::updateTransition(const HMM_PP::InferredDataFit* infer) {
	const uint n = infer->getData()->getNumObservations();
	if (!infer->hasDigamma())
		throw new Exception("Cannot updateTransition without having calculated digammas");

	for (uint i = 0; i < _params->getNumHiddenStates(); i++) {
		real denom = 0;
		for (uint t = 0; t < n-1; t++)
			denom += infer->getGamma(t,i);

		for (uint j = 0; j < _params->getNumHiddenStates(); j++) {
			real numer = 0;
			for (uint t = 0; t < n-1; t++)
				numer += infer->getDigamma(t,i,j);

			_transitionMat(i,j) += (numer / denom).getValue();
		}
	}
}

void HMM_PP::HomogeneousModelParams::scale(const BaseReal n) {
	_transitionMat /= n;
}

ostream& HMM_PP::HomogeneousModelParams::print(ostream& stream) const {
	stream << "A" << endl
			<< _transitionMat << endl;
	return stream;
}
