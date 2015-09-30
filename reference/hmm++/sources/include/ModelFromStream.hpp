#ifndef __MODEL_FROM_STREAM_H__
#define __MODEL_FROM_STREAM_H__

#include "Model.hpp"

#include <vector>
#include <utility>
using namespace std;

namespace HMM_PP {

	// The HMM
	class ModelFromStream : public Model {

	public:
		ModelFromStream(ModelParams* modelParams);
		virtual ~ModelFromStream();

		// load/clear (observed) data:
		virtual void loadAllData();

	protected:
		virtual uint addData(Data* d);

	public:
		// Get observed and (inferred hidden) sequence data from HMM:
		virtual const Data* getData(const uint i);
		virtual InferredDataFit* getDataFitResult(const uint i);

	protected:
		virtual void setDataFitResult(InferredDataFit* idf, const uint ind);

	public:
		virtual uint getNumIndividuals() { return _sequences.size(); }

	protected:
		// all attached sequences and the inference results for them:
		typedef pair<Data*, InferredDataFit*> DataAndFit;
		vector<DataAndFit> _sequences;
	};
}

#endif
