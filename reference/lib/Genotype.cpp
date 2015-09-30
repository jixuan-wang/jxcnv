#include "Genotype.hpp"
#include "SampleGenotypeQualsCalculator.hpp"
#include "AlleleQuals.hpp"

#include <Data.hpp>
#include <Model.hpp>

#include <sstream>
using namespace std;

const BaseReal XHMM::Genotype::MAX_PL = 255;

#define INIT \
		_sample(genotypeCalc->getSample()), \
		_mergedTargets(range == NULL ? NULL : new Interval(genotypeCalc->getMergedInterval(range->getStartTargIndex(), range->getStopTargIndex()))), \
		_targRange(range == NULL ? NULL : new TargetRange(*range)), \
		_meanRD(range == NULL ? NULL : genotypeCalc->calcMeanRD(range->getStartTargIndex(), range->getStopTargIndex())), \
		_meanOrigRD(meanOrigRD), \
		_nonDiploidScore(range == NULL ? NULL : new BaseReal(genotypeCalc->calcNonDiploidScore(range->getStartTargIndex(), range->getStopTargIndex()))), \
		_diploidScore(range == NULL ? NULL : new BaseReal(genotypeCalc->calcDiploidScore(range->getStartTargIndex(), range->getStopTargIndex())))

XHMM::Genotype::Genotype(const SampleGenotypeQualsCalculator* genotypeCalc, const TargetRange* range, const uint refType, const vector<uint>& nonRefTypes, const HMM_PP::ModelParams::DataVal* meanOrigRD, const RealThresh* callQualThresh)
: INIT {
	createNonRefAllelesAndCallGenotype(genotypeCalc, range, refType, nonRefTypes, callQualThresh);
}

XHMM::Genotype::Genotype(const SampleGenotypeQualsCalculator* genotypeCalc, const TargetRange* range, const uint refType, const uint nonRefType, const HMM_PP::ModelParams::DataVal* meanOrigRD, const RealThresh* callQualThresh)
: INIT {
	vector<uint> nonRefTypes;
	nonRefTypes.push_back(nonRefType);
	createNonRefAllelesAndCallGenotype(genotypeCalc, range, refType, nonRefTypes, callQualThresh);
}

XHMM::Genotype::~Genotype() {
	if (_mergedTargets != NULL)
		delete _mergedTargets;

	if (_targRange != NULL)
		delete _targRange;

	if (_meanRD != NULL)
		delete _meanRD;

	if (_meanOrigRD != NULL)
		delete _meanOrigRD;

	if (_nonDiploidScore != NULL)
		delete _nonDiploidScore;

	if (_diploidScore != NULL)
		delete _diploidScore;

	if (_nonRefAlleles != NULL) {
		for (map<uint, AlleleQuals*>::const_iterator it = _nonRefAlleles->begin(); it != _nonRefAlleles->end(); ++it)
			delete it->second;
		delete _nonRefAlleles;
	}

	if (_allAllelesToPLscores != NULL)
		delete _allAllelesToPLscores;
}

void XHMM::Genotype::createNonRefAllelesAndCallGenotype(const SampleGenotypeQualsCalculator* genotyper, const TargetRange* range, const uint refType, const vector<uint>& nonRefTypes, const RealThresh* callQualThresh) {
	_nonRefAlleles = NULL;
	_allAllelesToPLscores = NULL;

	if (range != NULL) {
		real refLikelihood = genotyper->calcLikelihoodGivenExactEvent(range->getStartTargIndex(), range->getStopTargIndex(), refType);

		map<uint, real>* allAllelesToLikelihoods = new map<uint, real>();
		(*allAllelesToLikelihoods)[refType] = refLikelihood;

		_nonRefAlleles = new map<uint, AlleleQuals*>();
		real maxLikelihoodGivenEvent = refLikelihood;
		for (vector<uint>::const_iterator typeIt = nonRefTypes.begin(); typeIt != nonRefTypes.end(); ++typeIt) {
			const uint nonRefType = *typeIt;

			AlleleQuals* aq = new AlleleQuals(genotyper, range->getStartTargIndex(), range->getStopTargIndex(), nonRefType);
			(*_nonRefAlleles)[nonRefType] = aq;

			const real& likelihoodGivenEvent = aq->getLikelihoodGivenEvent();
			(*allAllelesToLikelihoods)[nonRefType] = likelihoodGivenEvent;
			if (likelihoodGivenEvent > maxLikelihoodGivenEvent)
				maxLikelihoodGivenEvent = likelihoodGivenEvent;
		}

		_allAllelesToPLscores = new map<uint, BaseReal>();
		for (map<uint, real>::const_iterator likIt = allAllelesToLikelihoods->begin(); likIt != allAllelesToLikelihoods->end(); ++likIt) {
			real lik = likIt->second;
			real relLik = lik / maxLikelihoodGivenEvent;
			(*_allAllelesToPLscores)[likIt->first] = SampleGenotypeQualsCalculator::probToPhredScale(relLik, MAX_PL);
		}
		delete allAllelesToLikelihoods;
	}

	callGenotype(callQualThresh);
}

XHMM::Genotype::callType XHMM::Genotype::callGenotype(const RealThresh* callQualThresh) {
	_callType = NO_CALL;
	_calledAllele = NULL;

	if (callQualThresh != NULL && _nonRefAlleles != NULL && getDiploidScore() != NULL) {
		const AlleleQuals* maxSomeEventAllele = NULL;
		BaseReal max_SQ = 0;

		for (map<uint, AlleleQuals*>::const_iterator allIt = _nonRefAlleles->begin(); allIt != _nonRefAlleles->end(); ++allIt) {
			const AlleleQuals* all = allIt->second;
			BaseReal SQ_type = all->getHaveSomeEventScore();
			if (SQ_type >= max_SQ) {
				max_SQ = SQ_type;
				maxSomeEventAllele = all;
			}
		}

		if (callQualThresh->passScore(*getDiploidScore()) && callQualThresh->failScore(max_SQ)) // sample has the REF allele
			_callType = REFERENCE;
		else if (callQualThresh->passScore(maxSomeEventAllele->getHaveExactEventScore()) && callQualThresh->failScore(maxSomeEventAllele->getDontHaveAnyEventScore())) { // sample has argmaxType ALT allele
			_callType = NON_REFERENCE;
			_calledAllele = maxSomeEventAllele;
		}
	}

	return _callType;
}

const BaseReal* XHMM::Genotype::getPL(const uint type) const {
	if (_allAllelesToPLscores == NULL)
		return NULL;

	map<uint, BaseReal>::const_iterator findIt = _allAllelesToPLscores->find(type);
	if (findIt == _allAllelesToPLscores->end()) {
		stringstream str;
		str << "Cannot lookup PL for type " << type;
		throw new Exception(str.str());
	}

	return &(findIt->second);
}
