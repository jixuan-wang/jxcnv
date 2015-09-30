#include "CNVmodelParams.hpp"
#include "SampleGenotypeQualsCalculator.hpp"
#include "ReadDepthMatrixLoader.hpp"
#include "xhmmOutputManager.hpp"

#include <utils.hpp>
#include <Data.hpp>
#include <Model.hpp>
#include <Exception.hpp>
#include <Simplex.hpp>
#include <PreciseNonNegativeReal.hpp>
#include <InferredDataFit.hpp>

#include <cmath>
#include <sstream>
using namespace std;

const BaseReal XHMM::SampleGenotypeQualsCalculator::PHRED_SCALE_FACTOR = 10.0;

string XHMM::SampleGenotypeQualsCalculator::getSample() const {
	return getDataFit()->getData()->getId();
}

/*
 * CNV-specific HMM queries
 */
real XHMM::SampleGenotypeQualsCalculator::calcPosteriorRatioToDiploid(const uint t1, const uint t2, const uint type) const {
	real ratio = _sampleDataFit->calcPosteriorNumerator(t1, t2, type);
	ratio /= _sampleDataFit->calcPosteriorNumerator(t1, t2, CNVmodelParams::DIPLOID);
	return ratio;
}

pair<real, real> XHMM::SampleGenotypeQualsCalculator::calcCNVedgePosteriors(const uint t1, const uint t2, const uint type) const {
	real startProb;
	if (t1 > 0) {
		if (_sampleDataFit->hasDigamma()) // Use pre-calculated values:
			startProb = _sampleDataFit->getDigamma(t1 - 1, CNVmodelParams::DIPLOID, type);
		else
			startProb = _sampleDataFit->calcPosterior(t1 - 1, getStartBreakpointTransitionPath(type));
	}
	else
		startProb = 0;

	real stopProb;
	if (t2 < _sampleDataFit->getData()->getNumObservations() - 1) {
		if (_sampleDataFit->hasDigamma()) // Use pre-calculated values:
			stopProb = _sampleDataFit->getDigamma(t2, type, CNVmodelParams::DIPLOID);
		else
			stopProb = _sampleDataFit->calcPosterior(t2, getStopBreakpointTransitionPath(type));
	}
	else
		stopProb = 0;

	return pair<real, real>(startProb, stopProb);
}

pair<real, real> XHMM::SampleGenotypeQualsCalculator::calcCNVedgePosteriorComplements(const uint t1, const uint t2, const uint type) const {
	real startProbComplement;
	if (t1 > 0)
		startProbComplement = _sampleDataFit->calcPosteriorComplement(t1 - 1, getStartBreakpointTransitionPath(type));
	else
		startProbComplement = 1;

	real stopProbComplement;
	if (t2 < _sampleDataFit->getData()->getNumObservations() - 1)
		stopProbComplement = _sampleDataFit->calcPosteriorComplement(t2, getStopBreakpointTransitionPath(type));
	else
		stopProbComplement = 1;

	return pair<real, real>(startProbComplement, stopProbComplement);
}

vector<uint> XHMM::SampleGenotypeQualsCalculator::getStartBreakpointTransitionPath(const uint type) {
	vector<uint> startPath;

	startPath.push_back(CNVmodelParams::DIPLOID);
	startPath.push_back(type);

	return startPath;
}

vector<uint> XHMM::SampleGenotypeQualsCalculator::getStopBreakpointTransitionPath(const uint type) {
	vector<uint> stopPath;

	stopPath.push_back(type);
	stopPath.push_back(CNVmodelParams::DIPLOID);

	return stopPath;
}

set<uint>* XHMM::SampleGenotypeQualsCalculator::getOppositeStates(const uint type) {
	uint oppositeType;

	if (type == CNVmodelParams::DEL)
		oppositeType = CNVmodelParams::DUP;
	else if (type == CNVmodelParams::DUP)
		oppositeType = CNVmodelParams::DEL;
	else
		oppositeType = CNVmodelParams::DIPLOID;

	set<uint>* oppositeStates = new set<uint>();
	oppositeStates->insert(oppositeType);

	return oppositeStates;
}

real XHMM::SampleGenotypeQualsCalculator::calcPosteriorNumeratorCalledTypeNoOppositeType(const uint t1, const uint t2, const uint type) const {
	set<uint>* excludedStates = getOppositeStates(type);

	// Calculates the posterior probability: Pr(x \in (type, DIPLOID)+)
	real posteriorExcludedNumerator = _sampleDataFit->calcPosteriorNumeratorAssignmentsExcludingStates(t1, t2, *excludedStates);
	delete excludedStates;

	// Subtract out the all-DIPLOID call:
	posteriorExcludedNumerator -= _sampleDataFit->calcPosteriorNumerator(t1, t2, CNVmodelParams::DIPLOID);

	return posteriorExcludedNumerator;
}

real XHMM::SampleGenotypeQualsCalculator::calcPosteriorCalledTypeNoOppositeType(const uint t1, const uint t2, const uint type) const {
	real posterior = calcPosteriorNumeratorCalledTypeNoOppositeType(t1, t2, type);
	posterior /= _sampleDataFit->calcLikelihood();
	return posterior;
}

real XHMM::SampleGenotypeQualsCalculator::calcPosteriorComplementCalledTypeNoOppositeType(const uint t1, const uint t2, const uint type) const {
	real posteriorNumerator = calcPosteriorNumeratorCalledTypeNoOppositeType(t1, t2, type);
	real likelihood = _sampleDataFit->calcLikelihood();

	real posteriorComplement = likelihood - posteriorNumerator;
	posteriorComplement /= likelihood;

#if defined(DEBUG)
	set<uint>* excludedStates = getOppositeStates(type);
	real posteriorComplement_2 = _sampleDataFit->calcPosteriorComplementAssignmentsExcludingStates(t1, t2, *excludedStates);
	delete excludedStates;

	posteriorComplement_2 += _sampleDataFit->calcPosterior(t1, t2, CNVmodelParams::DIPLOID);

	if (abs(posteriorComplement.getLog10Value() - posteriorComplement_2.getLog10Value()) > 1e-1) {
		stringstream str;
		str << "calcPosteriorComplementCalledTypeNoOppositeType - precision error: " << -posteriorComplement.getLog10Value() << " != " << -posteriorComplement_2.getLog10Value();
		XHMM::xhmmOutputManager::printWarning(str.str());
	}
#endif

	return posteriorComplement;
}

/*
 * CNV event likelihoods
 */
real XHMM::SampleGenotypeQualsCalculator::calcLikelihoodGivenExactEvent(const uint t1, const uint t2, const uint state) const {
	return _sampleDataFit->calcLikelihoodGivenEvent(t1, t2, state);
}

/*
 * CNV event scores
 */
BaseReal XHMM::SampleGenotypeQualsCalculator::calcHaveExactEventScore(const uint t1, const uint t2, const uint type) const {
	return calcScoreFromErrorProb(_sampleDataFit->calcPosteriorComplement(t1, t2, type));
}

BaseReal XHMM::SampleGenotypeQualsCalculator::calcHaveSomeEventScore(const uint t1, const uint t2, const uint type) const {
	return calcScoreFromErrorProb(calcPosteriorComplementCalledTypeNoOppositeType(t1, t2, type));
}

BaseReal XHMM::SampleGenotypeQualsCalculator::calcDontHaveAnyEventScore(const uint t1, const uint t2, const uint type) const {
	set<uint> excludedStates;
	excludedStates.insert(type);

	return calcScoreFromErrorProb(_sampleDataFit->calcPosteriorComplementAssignmentsExcludingStates(t1, t2, excludedStates));
}

BaseReal XHMM::SampleGenotypeQualsCalculator::calcNonDiploidScore(const uint t1, const uint t2) const {
	return calcScoreFromErrorProb(_sampleDataFit->calcPosterior(t1, t2, CNVmodelParams::DIPLOID));
}

BaseReal XHMM::SampleGenotypeQualsCalculator::calcDiploidScore(const uint t1, const uint t2) const {
	return calcHaveExactEventScore(t1, t2, CNVmodelParams::DIPLOID);
}

pair<BaseReal, BaseReal> XHMM::SampleGenotypeQualsCalculator::calcStartStopScores(const uint t1, const uint t2, const uint type) const {
	pair<real, real> startStopProbComplements = calcCNVedgePosteriorComplements(t1, t2, type);

	return pair<BaseReal, BaseReal>(calcScoreFromErrorProb(startStopProbComplements.first), calcScoreFromErrorProb(startStopProbComplements.second));
}

BaseReal XHMM::SampleGenotypeQualsCalculator::calcRatioScore(const uint t1, const uint t2, const uint type) const {
	return - calcScoreFromErrorProb(calcPosteriorRatioToDiploid(t1, t2, type));
}

BaseReal XHMM::SampleGenotypeQualsCalculator::calcScoreFromErrorProb(const real& errorProb) const {
	return XHMM::SampleGenotypeQualsCalculator::probToPhredScale(errorProb, _maxScoreVal);
}

BaseReal XHMM::SampleGenotypeQualsCalculator::probToPhredScale(const real& prob, const BaseReal maxPhredScaleVal) {
	BaseReal score = - PHRED_SCALE_FACTOR * prob.getLog10Value();
	if (score <= 0 && score >= - real::REAL_EPSILON) // set rounding errors [and -0] to 0
		score = 0;

	if (maxPhredScaleVal > 0)
		HMM_PP::capAbsoluteValue(score, maxPhredScaleVal);

	return score;
}

HMM_PP::ModelParams::DataVal* XHMM::SampleGenotypeQualsCalculator::calcMeanRD(const uint t1, const uint t2) const {
	return _model->calcRepresentativeDataVal(_sampleDataFit->getData(), t1, t2);
}

XHMM::Interval XHMM::SampleGenotypeQualsCalculator::getMergedInterval(const uint t1, const uint t2) const {
	return _dataLoader->getTarget(t1) + _dataLoader->getTarget(t2);
}
