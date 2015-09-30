#ifndef __ALLELE_QUALS_H__
#define __ALLELE_QUALS_H__

#include <params.hpp>
#include <ModelParams.hpp>
#include <PreciseNonNegativeReal.hpp>

#include <utility>
using namespace std;

namespace XHMM {
	class SampleGenotypeQualsCalculator;

	class AlleleQuals {

	public:
		AlleleQuals(const SampleGenotypeQualsCalculator* genotypeCalc, const uint t1, const uint t2, const uint type);
		~AlleleQuals();

		uint getCNVtype() const { return _cnvType; }

		BaseReal getHaveExactEventScore() const { return _eventScore; }
		BaseReal getHaveSomeEventScore() const { return _someEventScore; }
		BaseReal getDontHaveAnyEventScore() const { return _noEventScore; }

		const pair<BaseReal, BaseReal>& getStartStopScores() const { return _startStopScores; }

		BaseReal getRatioScore() const { return _ratioScore; }

		const real& getLikelihoodGivenEvent() const { return _likelihoodGivenEvent; }

	private:
		uint _cnvType;

		BaseReal _eventScore;
		BaseReal _someEventScore;
		BaseReal _noEventScore;

		pair<BaseReal, BaseReal> _startStopScores;

		BaseReal _ratioScore;

		real _likelihoodGivenEvent;
	};

}

#endif
