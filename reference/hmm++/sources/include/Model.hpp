#ifndef __MODEL_H__
#define __MODEL_H__

#include "utils.hpp"
#include "Exception.hpp"
#include "NamedVector.hpp"
#include "NamedMatrix.hpp"
#include "PreciseNonNegativeReal.hpp"
#include "ModelParams.hpp"

#include <string>
#include <map>
#include <vector>
#include <istream>
#include <utility>
using namespace std;

namespace HMM_PP {
	class Data;
	class InferredDataFit;

	// The HMM
	class Model {

	public:
		Model(ModelParams* modelParams);
		virtual ~Model();

		inline ModelParams* getModelParams() { return _modelParams; }
		inline const ModelParams* getModelParams() const { return _modelParams; }

		inline Data* getNewData() { return _modelParams->getNewData(); }
		inline ModelParams::DataVal* calcRepresentativeDataVal(const HMM_PP::Data* d, const uint t1, const uint t2) const { return _modelParams->calcRepresentativeDataVal(d, t1, t2); }

		// load/clear (observed) data:
		virtual void loadAllData();

	protected:
		virtual uint addData(Data* d) = 0;

	public:
		// Get observed and (inferred hidden) sequence data from HMM:
		virtual const Data* getData(const uint i) = 0;
		virtual InferredDataFit* getDataFitResult(const uint i) = 0;

		// MUST be called after getData() to release Data:
		virtual void finishData(Data* d) const {}

		// MUST be called after getDataFitResult() to release finishDataAndFit():
		virtual void finishDataAndFit(InferredDataFit* idf) const {}

	protected:
		virtual void setDataFitResult(InferredDataFit* idf, const uint ind) = 0;

	private:
		void _internal_setDataFitResult(InferredDataFit* idf, const uint ind);

	protected:
		pair<InferredDataFit*, real> fitSequence(const Data* seq, bool calcDigammas, bool viterbi, bool timeEvents);

	public:
		// HMM operations
		virtual real fit(bool calcDigammas, bool viterbi, bool timeEvents = false);

		void baumWelch();
		bool nelderMead();

		void reestimateTransitionAndEmissionNoScaling();
		void reestimateInitProbs();

		// queries
		void displayFittedSequence(ostream& stream);

		inline uint getNumHiddenStates() const { return _modelParams->getNumHiddenStates(); }

		virtual uint getNumIndividuals() = 0;

		// Misc
		inline friend ostream& operator<<(ostream& stream, Model const& m) {
			return m.printModel(stream);
		}

		virtual ostream& printModel(ostream& stream) const;
		virtual ostream& printExplicitModel(ostream& stream, const uint maxT = 1) const;

		// transition probability function
		BaseReal transitionFunc(const uint i, const uint j, const uint t1, const uint t2) const;

		// convenience methods for transition probability function:
		HMM_PP::NamedMatrix<real>* transitionMatrix(const uint t1, const uint t2) const;
		HMM_PP::NamedMatrix<real>* transitionMatrixTranspose(const uint t1, const uint t2) const;

		// emission function
		real emissionFunc(const Data* seq, const uint hidden_state, const uint t) const;
		vector<real>* emissionFunc(const HMM_PP::Data* seq, const uint t) const;

	protected:
		// pointer to 'a' and 'b' specification
		ModelParams* _modelParams;
	};
}

#endif
