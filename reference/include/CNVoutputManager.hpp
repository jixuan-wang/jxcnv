#ifndef __CNV_OUTPUT_MANAGER_H__
#define __CNV_OUTPUT_MANAGER_H__

#include <params.hpp>
#include <DataOutputter.hpp>
#include <InferredDataFit.hpp>

#include <cstddef>
using namespace std;

namespace HMM_PP {
	class Model;
	class Data;
}

namespace XHMM {
	class ReadDepthMatrixLoader;

	class CNVoutputManager : public HMM_PP::DataOutputter<HMM_PP::InferredDataFit> {

	public:
		CNVoutputManager(ReadDepthMatrixLoader* origDataLoader = NULL);
		virtual ~CNVoutputManager() = 0;

		virtual void setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader);

	protected:
		HMM_PP::Model* _model;
		const ReadDepthMatrixLoader* _dataLoader;

		HMM_PP::Data* getOrigDataForNextDataFitOutput(const HMM_PP::InferredDataFit* idf);
		bool hasOrigData() const { return _origDataLoader != NULL; }

	private:
		ReadDepthMatrixLoader* _origDataLoader;
	};

}

#endif
