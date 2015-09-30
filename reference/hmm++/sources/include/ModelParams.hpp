#ifndef __MODEL_PARAMS_H__
#define __MODEL_PARAMS_H__

#include "NamedVector.hpp"
#include "NamedMatrix.hpp"
#include "utils.hpp"
#include "params.hpp"
#include "PreciseNonNegativeReal.hpp"

#include <string>
#include <sstream>
using namespace std;

namespace HMM_PP {
	class Data;
	class Model;
	class InferredDataFit;

	class ModelParams {

	protected:
		class HiddenStateParams;
		ModelParams(HiddenStateParams* params);

	public:
		virtual ~ModelParams() = 0;

		inline uint getNumHiddenStates() const { return _params->getNumHiddenStates(); }

		// Get/calculate values
		virtual BaseReal transitionFunc(const uint i, const uint j, const uint t1, const uint t2) const = 0;

		virtual real emissionFunc(const HMM_PP::Data* sequence, const uint i, const uint t) const = 0;
		virtual vector<real>* emissionFunc(const HMM_PP::Data* sequence, const uint t) const;

		// Specify 'a' and 'b' by a list of parameters (from simplex)
		virtual BaseRealVec getVariableParams() const = 0;
		virtual void setVariableParams(const BaseRealVec& v) = 0;

		// Updates
		virtual void updateTransition(const HMM_PP::InferredDataFit*) = 0;
		virtual void updateEmission(const HMM_PP::InferredDataFit*) = 0;
		virtual void clear() = 0;
		virtual void scale(const BaseReal n) = 0;

		virtual Data* getNewData() const = 0;
		virtual bool hasNextLoadData() const = 0;
		virtual Data* loadNextData() const = 0;

		const BaseRealVec& getInitProbs() const { return _params->_initProbs; }
		void setStartingProbsAndNormalize(const BaseRealVec& x) { _params->setStartingProbsAndNormalize(x); }

		const string& state(const uint i) const { return _params->state(i); }
		uint state(const string& s) const { return _params->state(s); }

		vector<uint> state(const vector<string>& s) const { return _params->state(s); }
		vector<string> state(const vector<uint>& s) const { return _params->state(s); }

		template<class RealType>
		void giveStateNames(NamedMatrix<RealType>& mat) const;

		template<class RealType>
		void giveStateNames(NamedVector<RealType>& vec) const;

		inline friend ostream& operator<<(ostream& stream, ModelParams const& mp) {
			return mp.print(stream);
		}

		class DataVal {
		public:
			DataVal() {}
			virtual ~DataVal() {}

			inline friend ostream& operator<<(ostream& stream, DataVal const& dv) {
				return dv.print(stream);
			}

			virtual ostream& print(ostream& stream) const = 0;
		};

		virtual DataVal* calcRepresentativeDataVal(const HMM_PP::Data* d, const uint t1, const uint t2) const = 0;

	protected:
		HiddenStateParams* _params;

		static real dnorm(const BaseReal x, const BaseReal u, const BaseReal v);
		static BaseReal univarNorm(const BaseReal x, const BaseReal u, const BaseReal s);

		class HiddenStateParams {
		public:
			HiddenStateParams() : _statemap(), _statermap(), _initProbs() {}
			HiddenStateParams(istream& stream);

			~HiddenStateParams() {}

			// name labels for hidden states:
			vector<string> _statemap;
			map<string, int> _statermap;

			// starting probs
			BaseRealVec _initProbs;

			void setNumHiddenStates(const uint k);
			void addState(const string& label);
			void relabelHiddenState(const uint i, const string& label);

			const string& state(const uint i) const { if (i >= _statemap.size()) { stringstream str; str << "Undefined state in model: " << i; throw new Exception(str.str()); } return _statemap[i]; }
			uint state(const string& s) const { if (_statermap.find(s) == _statermap.end()) throw new Exception("Undefined state in model: " + s); return _statermap.find(s)->second; }

			vector<uint> state(const vector<string>& s) const;
			vector<string> state(const vector<uint>& s) const;

			inline uint getNumHiddenStates() const { return _statemap.size(); }

			BaseRealVec getUniformStartProbs() const;
			void setStartingProbsAndNormalize(const BaseRealVec& x);
		};

		// display model
		virtual ostream& print(ostream& stream) const = 0;

		static const BaseReal LOG_10_E;
	};

	template<class RealType>
	void HMM_PP::ModelParams::giveStateNames(NamedMatrix<RealType>& mat) const {
		ullintPair dims = mat.getDims();
		if (dims.first != getNumHiddenStates() || dims.second != getNumHiddenStates())
			throw new Exception("Cannot assign names to incorrectly sized matrix");

		for (uint i = 0; i < getNumHiddenStates(); ++i) {
			const string& st = this->state(i);
			mat.setRowName(i, st);
			mat.setColName(i, st);
		}
	}

	template<class RealType>
	void HMM_PP::ModelParams::giveStateNames(NamedVector<RealType>& vec) const {
		if (vec.size() != getNumHiddenStates())
			throw new Exception("Cannot assign names to incorrectly sized vector");

		for (uint i = 0; i < getNumHiddenStates(); ++i) {
			const string& st = this->state(i);
			vec.setName(i, st);
		}
	}
}

#endif
