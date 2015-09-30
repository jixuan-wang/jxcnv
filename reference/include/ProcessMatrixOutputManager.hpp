#ifndef __PROCESS_MATRIX_OUTPUT_MANAGER_H__
#define __PROCESS_MATRIX_OUTPUT_MANAGER_H__

#include "xhmmCmdline.h"
#include "xhmmInputManager.hpp"

#include <params.hpp>

#include <string>
#include <list>
using namespace std;

namespace HMM_PP {
	template<class T> class NamedMatrix;
}

namespace XHMM {

	class ProcessMatrixOutputManager {

	public:
		ProcessMatrixOutputManager();
		~ProcessMatrixOutputManager();

		void processReadDepthMatrix(string rdFile, string outFile,
				list<string>* excludeTargetsFiles, list<string>* excludeChromosomeTargets, list<string>* excludeSamplesFiles,
				const uint minTargetSize, const uint maxTargetSize,
				const double minMeanTargetRD, const double maxMeanTargetRD,
				const double minSdTargetRD, const double maxSdTargetRD,
				const double minMeanSampleRD, const double maxMeanSampleRD,
				const double minSdSampleRD, const double maxSdSampleRD,
				bool scaleDataBySum, const enum_scaleDataBySumType& scaleDataBySumType, const double scaleDataBySumFactor,
				bool log10Data, const double pseudoCountValForLog10Input,
				bool centerData, const enum_centerType& centerType,
				bool zScoreData,
				string outputExcludedTargets, string outputExcludedSamples);

	private:
		HMM_PP::DoubleMat*
		filterTargetProperties(HMM_PP::DoubleMat* rdMat,
				const double minMeanTargetRD, const double maxMeanTargetRD,
				const double minSdTargetRD, const double maxSdTargetRD,
				xhmmInputManager::IntervalSet* excludedTargets);

		HMM_PP::DoubleMat*
		filterSampleProperties(HMM_PP::DoubleMat* rdMat,
				const double minMeanSampleRD, const double maxMeanSampleRD,
				const double minSdSampleRD, const double maxSdSampleRD,
				xhmmInputManager::StringSet* excludedSamples);
	};

}

#endif
