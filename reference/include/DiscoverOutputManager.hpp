#ifndef __DISCOVER_OUTPUT_MANAGER_H__
#define __DISCOVER_OUTPUT_MANAGER_H__

#include "xhmmInputManager.hpp"
#include "CNVoutputManager.hpp"

#include <params.hpp>
#include <PreciseNonNegativeReal.hpp>
#include <utils.hpp>

#include <utility>
#include <set>
#include <string>
#include <map>
using namespace std;

namespace HMM_PP {
	class Model;
}

namespace XHMM {
	class CNVmodelParams;
	class SampleGenotypeQualsCalculator;
	class AuxiliaryCNVoutputManager;
	class ReadDepthMatrixLoader;
	class PosteriorOutputManager;

	class DiscoverOutputManager : public CNVoutputManager {

	private:
		static const int DEFAULT_PRECISION = 2;

	public:
		DiscoverOutputManager
		(string outFile, AuxiliaryCNVoutputManager* auxOutput,
				string posteriorFile,
				BaseReal discoverSomeQualThreshold = 0, BaseReal maxScoreVal = real::REAL_INFINITY, int scoreDecimalPrecision = DEFAULT_PRECISION,
				ReadDepthMatrixLoader* origDataLoader = NULL);

		virtual ~DiscoverOutputManager();

		virtual void setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader);

		void printAllCNVs(HMM_PP::Model* model);

		virtual void printDataType(const HMM_PP::InferredDataFit* idf);

		typedef pair<BaseReal*, set<string>*> ThreshAndSamples;
		typedef map<Interval, ThreshAndSamples> IntervalToThreshMap;
		static IntervalToThreshMap* CNVtableToIntervals(xhmmInputManager::Table* cnvTable);

		static const string SAMPLE;
		static const string CNV;
		static const string INTERVAL;
		static const string EXACT_QUAL;

	private:
		HMM_PP::ostreamWriter* _outStream;
		AuxiliaryCNVoutputManager* _auxOutput;

		string _posteriorFile;
		PosteriorOutputManager* _posteriorOutput;

		const HMM_PP::PrecisionThreshold<BaseReal> _discoverSomeQualThreshold;
		BaseReal _maxScoreVal;
		int _scoreDecimalPrecision;

		set<uint>* _excludeSegmentTypes;

		void printCNVheader();
		void printCNV(SampleGenotypeQualsCalculator* genotyper, const uint t1, const uint t2, const uint type, HMM_PP::Data* origData);
	};

}

#endif
