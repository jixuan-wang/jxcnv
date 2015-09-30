#ifndef __GENOTYPE_H__
#define __GENOTYPE_H__

#include "Interval.hpp"
#include "TargetRange.hpp"

#include <params.hpp>
#include <ModelParams.hpp>

#include <utility>
#include <map>
#include <vector>
using namespace std;

namespace XHMM {
	class SampleGenotypeQualsCalculator;
	class AlleleQuals;

	class Genotype {

	public:
		typedef HMM_PP::PrecisionThreshold<BaseReal> RealThresh;

		Genotype(const SampleGenotypeQualsCalculator* genotypeCalc, const TargetRange* range, const uint refType, const vector<uint>& nonRefTypes, const HMM_PP::ModelParams::DataVal* meanOrigRD = NULL, const RealThresh* callQualThresh = NULL);
		Genotype(const SampleGenotypeQualsCalculator* genotypeCalc, const TargetRange* range, const uint refType, const uint nonRefType, const HMM_PP::ModelParams::DataVal* meanOrigRD = NULL, const RealThresh* callQualThresh = NULL);

		~Genotype();

		const string& getSample() const { return _sample; }

		const Interval* getMergedTargets() const { return _mergedTargets; }
		const TargetRange* getTargetRange() const { return _targRange; }

		const HMM_PP::ModelParams::DataVal* getMeanRD() const { return _meanRD; }
		const HMM_PP::ModelParams::DataVal* getMeanOrigRD() const { return _meanOrigRD; }

		const BaseReal* getNonDiploidScore() const { return _nonDiploidScore; }
		const BaseReal* getDiploidScore() const { return _diploidScore; }

		uint getNumNonRefAlleles() const { if (_nonRefAlleles == NULL) return 0; else return _nonRefAlleles->size(); }
		const AlleleQuals* getNonRefAllele(uint type) const { if (_nonRefAlleles == NULL) return NULL; else return (*_nonRefAlleles)[type]; }
		const map<uint, AlleleQuals*>* getNonRefAlleles() const { return _nonRefAlleles; }

		const BaseReal* getPL(const uint type) const;

		enum callType {
			NO_CALL,
			REFERENCE,
			NON_REFERENCE
		};

		callType callGenotype(const RealThresh* callQualThresh);

		bool isNoCall() const { return _callType == NO_CALL; }
		bool isReference() const { return _callType == REFERENCE; }
		bool isNonReference() const { return _callType == NON_REFERENCE; }

		const AlleleQuals* getNonRefCall() const { return _calledAllele; }

		static const BaseReal MAX_PL;

	private:
		string _sample;
		Interval* _mergedTargets;
		TargetRange* _targRange;

		const HMM_PP::ModelParams::DataVal* _meanRD;
		const HMM_PP::ModelParams::DataVal* _meanOrigRD;

		BaseReal* _nonDiploidScore;
		BaseReal* _diploidScore;

		map<uint, AlleleQuals*>* _nonRefAlleles;
		map<uint, BaseReal>* _allAllelesToPLscores;

		callType _callType;
		const AlleleQuals* _calledAllele;

		void createNonRefAllelesAndCallGenotype(const SampleGenotypeQualsCalculator* genotypeCalc, const TargetRange* range, const uint refType, const vector<uint>& nonRefTypes, const RealThresh* callQualThresh);
	};

}

#endif
