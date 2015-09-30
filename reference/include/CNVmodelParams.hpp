#ifndef __CNV_MODEL_PARAMS_H__
#define __CNV_MODEL_PARAMS_H__

#include <utils.hpp>
#include <UnivariateQuantitativeModelParams.hpp>
#include <params.hpp>
#include <NamedMatrix.hpp>
#include <PreciseNonNegativeReal.hpp>
#include <UnivariateQuantitativeData.hpp>

#include <utility>
#include <vector>
using namespace std;

namespace HMM_PP {
	class Data;
}

namespace XHMM {
	class ReadDepthMatrixLoader;

	class CNVmodelParams : public HMM_PP::UnivariateQuantitativeModelParams {

	public:
		CNVmodelParams(istream& stream, vector<uint>* interTargetDistances = new vector<uint>());
		CNVmodelParams(istream& stream, ReadDepthMatrixLoader* dataLoader);

		virtual ~CNVmodelParams();

	private:
		void init(istream& stream);

	protected:
		static vector<uint>* calcInterTargetDistances(const ReadDepthMatrixLoader* dataLoader);

		void constructFixedParams();
		void constructVariableParamsFromStream(istream& stream);

	public:
		// Updates
		virtual HMM_PP::BaseRealVec getVariableParams() const;
		virtual void setVariableParams(const HMM_PP::BaseRealVec& v);

		// Get/calculate values
		virtual BaseReal transitionFunc(const uint i, const uint j, const uint t1, const uint t2) const;
	protected:
		BaseReal transitionProb(const uint i, const uint j, const uint distInBases) const;
	public:
		HMM_PP::NamedMatrix<real>* transitionMatrix(const uint distInBases, bool assignNames = true) const;

		/*
		 * If not using Baum-Welch, the following four functions do not actually need to be implemented
		 */
		virtual void clear() { throw new Exception("should not be used"); }
		virtual void updateTransition(const HMM_PP::InferredDataFit*) { throw new Exception("should not be used"); }
		virtual void updateEmission(const HMM_PP::InferredDataFit*) { throw new Exception("should not be used"); }
		virtual void scale(const BaseReal n) { throw new Exception("should not be used"); }

	public:
		enum CNVtype {
			DEL = 0,
			DIPLOID,
			DUP,
			NUM_CNV_TYPES
		};

		static const vector<uint> NON_DIPLOID_TYPES;

		// All types, with DIPLOID first:
		static const vector<uint> ALL_CN_TYPES;

	private:
		static vector<uint> getNonDiploidTypes();

		// Get all types, with DIPLOID first:
		static vector<uint> getAllTypes();

	private:
		// auxiliary data: consecutive inter-target distances
		vector<uint>* _interTargetDistances;

		// fixed: 3 hidden states

		// transition matrix parameterized by 3 variables + 1 data-dependent variable (i.e. inhomogeneous matrix)
		// emission matrix is free (3 means, 3 variances)

		// 'a' base transition matrix (augmented in non-homogeneous manner in transitionFunc())
		BaseReal _p;      // Probability of starting a del CNV = Prob. of starting a dup CNV
		BaseReal _e;      // Expected CNV length (# of targets)
		BaseReal _d;      // Expected inter-target distance within CNVs
		BaseReal _lambda; // Parameter for exponential decay of inter-target distance distribution [ == 1 / _d ]

		// _p, _e, _d:
		static const uint NUM_TRANS_MAT_PARAMS_TO_READ = 3;

		inline const HMM_PP::BaseRealVec& getDiploidTransProbs() const {
			return _transitionMat[DIPLOID];
		}

	protected:
		// display model
		virtual ostream& print(ostream& stream) const;
	};
}

#endif
