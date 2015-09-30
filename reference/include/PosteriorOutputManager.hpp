#ifndef __POSTERIOR_OUTPUT_MANAGER_H__
#define __POSTERIOR_OUTPUT_MANAGER_H__

#include <params.hpp>
#include <utils.hpp>

#include <utility>
#include <vector>
using namespace std;

namespace HMM_PP {
	class Model;
}

namespace XHMM {
	class ReadDepthMatrixLoader;

	class PosteriorOutputManager {

	public:
		PosteriorOutputManager(string posteriorFile, const vector<uint>& types, const HMM_PP::Model* model, const ReadDepthMatrixLoader* readDepthLoader);

		~PosteriorOutputManager();

		void printSamplePosteriors(const HMM_PP::InferredDataFit* idf);

	private:
		typedef pair<const uint, HMM_PP::ostreamWriter*> CNtypeAndWriter;
		list<CNtypeAndWriter> _CNtypes;
	};

}

#endif
