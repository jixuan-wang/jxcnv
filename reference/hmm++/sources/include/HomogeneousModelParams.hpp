#ifndef __HOMOGENEOUS_MODEL_PARAMS_H__
#define __HOMOGENEOUS_MODEL_PARAMS_H__

#include "utils.hpp"
#include "ModelParams.hpp"
#include "NamedMatrix.hpp"
#include "params.hpp"

namespace HMM_PP {

	class HomogeneousModelParams : virtual public ModelParams {

	public:
		HomogeneousModelParams(HiddenStateParams* params);
		virtual ~HomogeneousModelParams() = 0;

	protected:
		BaseRealVec getHomogeneousStationaryProbs() const;

		static void readTransMatParams(istream& stream, const uint numHiddenStates, BaseRealVec& paramValues);
		void getTransMatParams(BaseRealVec& v) const;
		void setTransMatFromParams(const BaseRealVec& v, uint& vInd);

	public:
		// Get/calculate values
		virtual BaseReal transitionFunc(const uint i, const uint j, const uint t1, const uint t2) const;

		virtual void clear();
		virtual void updateTransition(const HMM_PP::InferredDataFit*);
		virtual void scale(const BaseReal n);

	protected:
		BaseRealMat _transitionMat;

		// display model
		virtual ostream& print(ostream& stream) const;
	};
}

#endif
