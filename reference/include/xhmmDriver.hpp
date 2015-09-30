#ifndef __XHMM_DRIVER_H__
#define __XHMM_DRIVER_H__

#include "xhmmCmdline.h"

#include <string>
#include <istream>
#include <set>
using namespace std;

namespace HMM_PP {
	class Model;
	template<class DataType> class DataOutputter;
	class InferredDataFit;
}

namespace XHMM {
	class ReadDepthMatrixLoader;
	class DiscoverOutputManager;
	class AuxiliaryCNVoutputManager;

	class xhmmDriver {

	public:
		xhmmDriver(const gengetopt_args_info& args_info);
		~xhmmDriver();

		bool run();

	private:
		const gengetopt_args_info& _args_info;

		istream* _paramStream;
		string* _dbFile;

		HMM_PP::Model* _model;
		ReadDepthMatrixLoader* _dataLoader;

		ReadDepthMatrixLoader* _origDataLoader;
		AuxiliaryCNVoutputManager* _auxOutput;

		set<string>* _keepIDs;

		void runCommonInits();

		void setChrOrderFromReference();
		void readSampleIDsToKeep();

		// Command-line modes:
		bool prepareTargets();

		bool mergeGATKdepths();

		bool processMatrix();

		bool PCA();
		bool normalize();

		bool discover();
		bool genotype();
		bool mergeVCFs();

		bool printHMM();
		bool transition();

		bool createDB();

		istream* openParamFile();
		string* getDBfile();

		void createModelAndLoader(HMM_PP::DataOutputter<HMM_PP::InferredDataFit>* dataOuputter = NULL);
		void createOrigDataLoaderAndAuxOutput();
	};

}

#endif
