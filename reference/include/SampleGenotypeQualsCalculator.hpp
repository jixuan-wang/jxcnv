#ifndef __SAMPLE_GENOTYPE_QUALS_CALCULATOR_H__
#define __SAMPLE_GENOTYPE_QUALS_CALCULATOR_H__

#include "Interval.hpp"

#include <PreciseNonNegativeReal.hpp>
#include <ModelParams.hpp>

#include <set>
using namespace std;

namespace HMM_PP {
	class Data;
	class Model;
}

namespace XHMM {
	class ReadDepthMatrixLoader;

	class SampleGenotypeQualsCalculator {

	public:
		SampleGenotypeQualsCalculator(const HMM_PP::InferredDataFit* sampleDataFit, const HMM_PP::Model* model, BaseReal maxScoreVal, const ReadDepthMatrixLoader* dataLoader)
		: _sampleDataFit(sampleDataFit), _model(model), _maxScoreVal(maxScoreVal), _dataLoader(dataLoader) {}

		~SampleGenotypeQualsCalculator() {}

		const HMM_PP::InferredDataFit* getDataFit() const { return _sampleDataFit; }
		string getSample() const;

		/*
		 * CNV-specific HMM queries
		 */
		real calcPosteriorRatioToDiploid(const uint t1, const uint t2, const uint type) const;
		pair<real, real> calcCNVedgePosteriors(const uint t1, const uint t2, const uint type) const;
		pair<real, real> calcCNVedgePosteriorComplements(const uint t1, const uint t2, const uint type) const;

		real calcPosteriorNumeratorCalledTypeNoOppositeType(const uint t1, const uint t2, const uint type) const;
		real calcPosteriorCalledTypeNoOppositeType(const uint t1, const uint t2, const uint type) const;
		real calcPosteriorComplementCalledTypeNoOppositeType(const uint t1, const uint t2, const uint type) const;

		/*
		 * CNV event likelihoods
		 */
		real calcLikelihoodGivenExactEvent(const uint t1, const uint t2, const uint state) const;

		/*
		 * CNV event scores
		 */
		BaseReal calcHaveExactEventScore(const uint t1, const uint t2, const uint type) const;
		BaseReal calcHaveSomeEventScore(const uint t1, const uint t2, const uint type) const;
		BaseReal calcDontHaveAnyEventScore(const uint t1, const uint t2, const uint type) const;

		BaseReal calcNonDiploidScore(const uint t1, const uint t2) const;
		BaseReal calcDiploidScore(const uint t1, const uint t2) const;

		pair<BaseReal, BaseReal> calcStartStopScores(const uint t1, const uint t2, const uint type) const;

		BaseReal calcRatioScore(const uint t1, const uint t2, const uint type) const;

	private:
		const HMM_PP::InferredDataFit* _sampleDataFit;
		const HMM_PP::Model* _model;
		const BaseReal _maxScoreVal;
		const ReadDepthMatrixLoader* _dataLoader;

		static const BaseReal PHRED_SCALE_FACTOR;

		static vector<uint> getStartBreakpointTransitionPath(const uint type);
		static vector<uint> getStopBreakpointTransitionPath(const uint type);

		static set<uint>* getOppositeStates(const uint type);

		BaseReal calcScoreFromErrorProb(const real& errorProb) const;

	public:
		static BaseReal probToPhredScale(const real& prob, const BaseReal maxPhredScaleVal = 0);

		/*
		 * Extract other info about CNV
		 */
		HMM_PP::ModelParams::DataVal* calcMeanRD(const uint t1, const uint t2) const;

		Interval getMergedInterval(const uint t1, const uint t2) const;
	};

}

#endif
