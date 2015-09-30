#ifndef __MODEL_FROM_STREAM_NO_STORAGE_H__
#define __MODEL_FROM_STREAM_NO_STORAGE_H__

#include "Model.hpp"

namespace HMM_PP {
	class Data;
	class ModelParams;
	class InferredDataFit;
	template<class DataType> class DataOutputter;

	// The HMM
	class ModelFromStreamNoStorage : public Model {

	public:
		ModelFromStreamNoStorage(ModelParams* modelParams, DataOutputter<InferredDataFit>* dataOuputter);
		virtual ~ModelFromStreamNoStorage();

		// load/clear (observed) data:
		virtual void loadAllData() {}

	protected:
		virtual uint addData(Data* d) { return 0; }

	public:
		// Since the individuals CANNOT be retrieved from memory:
		virtual const Data* getData(const uint i) { return NULL; }
		virtual InferredDataFit* getDataFitResult(const uint i) { return NULL; }

	protected:
		virtual void setDataFitResult(InferredDataFit* idf, const uint ind) {}

	public:
		// Since the individuals CANNOT be retrieved from memory:
		virtual uint getNumIndividuals() { return 0; }

		virtual real fit(bool calcDigammas, bool viterbi, bool timeEvents = false);

	private:
		DataOutputter<InferredDataFit>* _dataOuputter;
	};
}

#endif
