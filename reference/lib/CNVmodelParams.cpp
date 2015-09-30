#include "CNVmodelParams.hpp"
#include "ReadDepthMatrixLoader.hpp"

#include <utils.hpp>
#include <Data.hpp>
#include <Model.hpp>
#include <Exception.hpp>
#include <Simplex.hpp>
#include <PreciseNonNegativeReal.hpp>

#include <cmath>
#include <sstream>
#include <list>
using namespace std;

const vector<uint> XHMM::CNVmodelParams::NON_DIPLOID_TYPES = CNVmodelParams::getNonDiploidTypes();
const vector<uint> XHMM::CNVmodelParams::ALL_CN_TYPES = CNVmodelParams::getAllTypes();

#define HEADER_SEPARATOR "***********************************************************************"

XHMM::CNVmodelParams::CNVmodelParams(istream& stream, vector<uint>* interTargetDistances)
: HMM_PP::ModelParams(new HiddenStateParams()),
  HMM_PP::UnivariateQuantitativeModelParams(HMM_PP::ModelParams::_params),
  _interTargetDistances(interTargetDistances), _p(0), _e(0), _d(0), _lambda(0) {
	init(stream);
}

XHMM::CNVmodelParams::CNVmodelParams(istream& stream, ReadDepthMatrixLoader* dataLoader)
: HMM_PP::ModelParams(new HiddenStateParams()),
  HMM_PP::UnivariateQuantitativeModelParams(HMM_PP::ModelParams::_params, dataLoader),
  _interTargetDistances(calcInterTargetDistances(dataLoader)), _p(0), _e(0), _d(0), _lambda(0) {
	init(stream);
}

void XHMM::CNVmodelParams::init(istream& stream) {
	constructFixedParams();
	constructVariableParamsFromStream(stream);
	_params->setStartingProbsAndNormalize(this->getDiploidTransProbs());

#if defined(DEBUG)
	giveStateNames(_transitionMat);
#endif

	cerr << *this << endl;
}

XHMM::CNVmodelParams::~CNVmodelParams() {
	delete _interTargetDistances;
}

vector<uint>* XHMM::CNVmodelParams::calcInterTargetDistances(const ReadDepthMatrixLoader* dataLoader) {
	vector<uint>* interTargetDistances = new vector<uint>(dataLoader->getNumTargets() - 1);

	for (uint t = 0; t < dataLoader->getNumTargets() - 1; ++t)
		(*interTargetDistances)[t] = dataLoader->getTarget(t+1).distance(dataLoader->getTarget(t));

	return interTargetDistances;
}

void XHMM::CNVmodelParams::constructFixedParams() {
	_params->setNumHiddenStates(NUM_CNV_TYPES);

	_params->relabelHiddenState(DEL,     "DEL");
	_params->relabelHiddenState(DIPLOID, "DIP");
	_params->relabelHiddenState(DUP,     "DUP");
}

void XHMM::CNVmodelParams::constructVariableParamsFromStream(istream& stream) {
	HMM_PP::BaseRealVec paramValues;

	for (uint i = 0; i < NUM_TRANS_MAT_PARAMS_TO_READ; i++) {
		BaseReal x;
		stream >> x;
		if (!stream)
			throw new Exception("Failed to read model transition matrix parameters from stream");
		paramValues.addElement(x);
	}

	UnivariateQuantitativeModelParams::readGaussianParams(stream, _params->getNumHiddenStates(), paramValues);

	setVariableParams(paramValues);
}

HMM_PP::BaseRealVec XHMM::CNVmodelParams::getVariableParams() const {
	HMM_PP::BaseRealVec v;

	v.addElement(_p);
	v.addElement(_e);
	v.addElement(_d / Interval::KB);

	getGaussianParams(v);

	return v;
}

void XHMM::CNVmodelParams::setVariableParams(const HMM_PP::BaseRealVec& v) {
	// expecting 3 parameters to determine the transition matrix
	if (v.size() < NUM_TRANS_MAT_PARAMS_TO_READ)
		throw new Exception("Internal error in CNVmodelParams::setVariableParams()");

	uint vInd = 0;
	_p = v[vInd++];
	_e = v[vInd++];
	_d = v[vInd++] * Interval::KB; // kb --> bp

	BaseReal q = 1 / _e;
	_lambda = 1 / _d; // stored, as accessed by f_trans()

	// Model parameters for transition matrix:

	//    p = probability of starting a del/dup
	//    e = expected length (mean number of TARGETS in each del/dup)
	//    d = mean distance (in bp) between targets in del/dup event

	// data-dependent factors
	//    distances = vector of distances (in bp) between consecutive pair of targets

	// Starting values :  p = 10^-8;
	//                    e = 6;
	//                    d = 70000;  // mean inter-target distance

	// 'base' matrix for 'a'
	// Will attenuate DEL and DUP states at i by distance between targets i and j, to revert
	// to DIPLOID state distribution, following exponential factor 'f'.  Done on a target-target
	// specific basis (non-homogeneous HMM)
	_transitionMat.setDims(_params->getNumHiddenStates(), _params->getNumHiddenStates(), 1/(BaseReal)_params->getNumHiddenStates());

	_transitionMat(DEL,DEL) = 1 - q;   _transitionMat(DEL,DIPLOID) = q;              _transitionMat(DEL,DUP) = 0;
	_transitionMat(DIPLOID,DEL) = _p;  _transitionMat(DIPLOID,DIPLOID) = 1 - 2 *_p;  _transitionMat(DIPLOID,DUP) = _p;
	_transitionMat(DUP,DEL) = 0;       _transitionMat(DUP,DIPLOID) = q;              _transitionMat(DUP,DUP) = 1 - q;

	setGaussianFromParams(v, vInd);
}

BaseReal XHMM::CNVmodelParams::transitionFunc(const uint i, const uint j, const uint t1, const uint t2) const {
	if (t2 != t1 + 1 || t1 >= _interTargetDistances->size()) {
		stringstream str;
		str << "XHMM::CNVmodelParams::transitionFunc: " << t2 << " = t2 != t1 + 1 = " << (t1 + 1) << ", or t1= " << t1 << " is out of bounds";
		throw new Exception(str.str());
	}
	// For a non-homogeneous transition matrix, calculate an interval-specific (t1->t2) transition: look up the distance in _targetDistances.
	uint distInBases = (*_interTargetDistances)[t1];
	return transitionProb(i, j, distInBases);
}

BaseReal XHMM::CNVmodelParams::transitionProb(const uint i, const uint j, const uint distInBases) const {
	/* Return a weighted mixture of the original and
	 * the 'infinite distance' transition (i.e., same as initial probability distribution)
	 */
	BaseReal f; // f == 1 - pexp(dist, lambda)
	if (distInBases == UINT_INFINITY)
		f = 0;
	else
		f = exp(-(distInBases * _lambda));

	return f * _transitionMat(i,j) + (1 - f) * getDiploidTransProbs()[j];
}

HMM_PP::NamedMatrix<real>* XHMM::CNVmodelParams::transitionMatrix(const uint distInBases, bool assignNames) const {
	HMM_PP::NamedMatrix<real>* a = new HMM_PP::NamedMatrix<real>();
	a->setDims(_params->getNumHiddenStates(), _params->getNumHiddenStates());

	for (uint i = 0; i < _params->getNumHiddenStates(); ++i) {
		for (uint j = 0; j < _params->getNumHiddenStates(); j++) {
			(*a)(i,j) = this->transitionProb(i, j, distInBases);
		}
	}

	if (assignNames)
		giveStateNames(*a);

	return a;
}

ostream& XHMM::CNVmodelParams::print(ostream& stream) const {
	stream
	<< endl
	<< HEADER_SEPARATOR << endl
	<< "Input CNV parameters file:" << endl
	<< HEADER_SEPARATOR << endl
	<< getVariableParams() << endl
	<< HEADER_SEPARATOR << endl << endl;

	stream
	<< HEADER_SEPARATOR << endl
	<< "xhmm parameters:" << endl
	<< HEADER_SEPARATOR << endl
	<< "Pr(start DEL) = Pr(start DUP) = " << _p << endl
	<< "Mean number of targets in CNV [geometric distribution] = " << _e << endl
	<< "Mean distance between targets within CNV [exponential decay] = " << (_d / Interval::KB) << " KB" << endl
	<< endl;

	for (uint type = DEL; type != NUM_CNV_TYPES; ++type)
		stream << state(type) << " read depth distribution ~ N(mean=" << _mean[type] << ", var=" << _var[type] << ")" << endl;

#if defined(DEBUG)
	cerr << endl << "*Base* transition matrix:" << endl << _transitionMat << endl;
	cerr << "[!UNUSED!] stationary distribution of *base* transition matrix:" << endl << getHomogeneousStationaryProbs() << endl;
#endif

	stream
	<< HEADER_SEPARATOR << endl;

	return stream;
}

vector<uint> XHMM::CNVmodelParams::getNonDiploidTypes() {
	vector<uint> types;

	for (uint type = 0; type < NUM_CNV_TYPES; ++type) {
		if (type != DIPLOID)
			types.push_back(type);
	}

	return types;
}

// Get all types, with DIPLOID first:
vector<uint> XHMM::CNVmodelParams::getAllTypes() {
	vector<uint> nonDiploidTypes = getNonDiploidTypes();
	list<uint> typesList(nonDiploidTypes.begin(), nonDiploidTypes.end());

	typesList.push_front(DIPLOID);

	return vector<uint>(typesList.begin(), typesList.end());
}
