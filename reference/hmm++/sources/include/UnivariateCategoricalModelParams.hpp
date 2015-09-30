#ifndef __UNIVARIATE_CATEGORICAL_MODEL_PARAMS_H__
#define __UNIVARIATE_CATEGORICAL_MODEL_PARAMS_H__

#include "utils.hpp"
#include "HomogeneousModelParams.hpp"
#include "NamedMatrix.hpp"
#include "params.hpp"
#include "UnivariateCategoricalData.hpp"
#include "ModelParamsData.hpp"

namespace HMM_PP {

	class UnivariateCategoricalModelParams : public HomogeneousModelParams, public ModelParamsData<UnivariateCategoricalData> {

	public:
		UnivariateCategoricalModelParams(istream& stream, HMM_PP::DataLoader<UnivariateCategoricalData>* dataLoader);
		virtual ~UnivariateCategoricalModelParams() {}

	protected:
		void constructFixedParamsFromStream(istream& stream);
		void constructVariableParamsFromStream(istream& stream);

	public:
		// Updates
		virtual BaseRealVec getVariableParams() const;
		virtual void setVariableParams(const BaseRealVec& v);

		// Get/calculate values
		virtual real emissionFunc(const HMM_PP::Data* sequence, const uint i, const uint t) const;
		virtual vector<real>* emissionFunc(const HMM_PP::Data* sequence, const uint t) const;

		virtual void clear();
		virtual void updateEmission(const HMM_PP::InferredDataFit*);
		virtual void scale(const BaseReal n);

	private:
		uint _numObservedStates;

		RealMat _emissionMat;

	protected:
		// display model
		virtual ostream& print(ostream& stream) const;

	public:
		class CategoricalDataVal : public DataVal {
		public:
			CategoricalDataVal(uint val) : _val(val) {}
			virtual ~CategoricalDataVal() {}

			virtual ostream& print(ostream& stream) const {
				stream << DataType::typeToCategory(_val);
				return stream;
			}

		private:
			uint _val;
		};

		virtual CategoricalDataVal* calcRepresentativeDataVal(const HMM_PP::Data* d, const uint t1, const uint t2) const;
	};
}

#endif
