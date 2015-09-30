#include "AlleleQuals.hpp"
#include "SampleGenotypeQualsCalculator.hpp"

#include <Data.hpp>
#include <Model.hpp>

#include <sstream>
using namespace std;

XHMM::AlleleQuals::AlleleQuals(const SampleGenotypeQualsCalculator* genotypeCalc, const uint t1, const uint t2, const uint type)
: _cnvType(type),
  _eventScore(genotypeCalc->calcHaveExactEventScore(t1, t2, type)),
  _someEventScore(genotypeCalc->calcHaveSomeEventScore(t1, t2, type)),
  _noEventScore(genotypeCalc->calcDontHaveAnyEventScore(t1, t2, type)),
  _startStopScores(genotypeCalc->calcStartStopScores(t1, t2, type)),
  _ratioScore(genotypeCalc->calcRatioScore(t1, t2, type)),
  _likelihoodGivenEvent(genotypeCalc->calcLikelihoodGivenExactEvent(t1, t2, type)) {
}

XHMM::AlleleQuals::~AlleleQuals() {
}
