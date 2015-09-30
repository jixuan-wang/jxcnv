#ifndef __GENOTYPE_OUTPUT_MANAGER_H__
#define __GENOTYPE_OUTPUT_MANAGER_H__

#include "ReadDepthMatrixLoader.hpp"
#include "CNVoutputManager.hpp"
#include "DiscoverOutputManager.hpp"
#include "Genotype.hpp"
#include "OnDiskMatrixTransposer.hpp"

#include <params.hpp>
#include <PreciseNonNegativeReal.hpp>

#include <ostream>
#include <vector>
#include <set>
#include <map>
#include <list>
using namespace std;

namespace HMM_PP {
	class Model;
}

namespace XHMM {
	class Interval;
	class CopyNumberVariant;
	class AuxiliaryCNVoutputManager;

	class GenotypeOutputManager : public CNVoutputManager {

	private:
		static const int DEFAULT_PRECISION = 2;

	public:
		GenotypeOutputManager(string outFile, const DiscoverOutputManager::IntervalToThreshMap* intervals, AuxiliaryCNVoutputManager* auxOutput, const BaseReal genotypeQualThresholdWhenNoExact, const bool genotypeSubsegments, const uint maxTargetsInSubsegment, ReadDepthMatrixLoader* origDataLoader = NULL, BaseReal maxScoreVal = real::REAL_INFINITY, int scoreDecimalPrecision = DEFAULT_PRECISION);
		virtual ~GenotypeOutputManager();

		virtual void setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader);

		virtual void printDataType(const HMM_PP::InferredDataFit* idf);

		static const string MISSING_VAL;

		static const string HEADER_ROW_PREFIX;

        static const string VCF_COLUMNS_ARRAY[];
        static list<string> VCF_COLUMNS;

        enum VCF_COLUMN_INDS {
        	CHROM = 0,
        	POS,
        	ID,
        	REF,
        	ALT,
        	QUAL,
        	FILTER,
        	INFO,
        	FORMAT
        };

	private:
		HMM_PP::ostreamWriter* _outStream;

		typedef pair<Genotype::RealThresh, set<string>*> ThreshAndSamples;
		typedef map<Interval, ThreshAndSamples> IntervalToPrecisionThreshMap;
		IntervalToPrecisionThreshMap* _intervals;

		AuxiliaryCNVoutputManager* _auxOutput;

		const BaseReal _genotypeQualThresholdWhenNoExact;

		const bool _genotypeSubsegments;
		const uint _maxTargetsInSubsegment;

		BaseReal _maxScoreVal;
		int _scoreDecimalPrecision;

		vector<string>* _samples;
		map<Interval, CopyNumberVariant*>* _allGenotypes;

		void printVCFheader();

	public:
		static void initOstream(ostream& stream);
		static void printSampleHeaderLine(ostream& outStream, const vector<string>* samples);

	private:
		void printVCFgenotypesAllIntervals();

		string generateGenotypeString(const Genotype* gt, const CopyNumberVariant* cnv) const;

		/*
		 * Transpose functions:
		 */
		string _vcfTransposerFile;
		OnDiskMatrixTransposer* _vcfTransposer;

		typedef map<string, string> SampleToGenotypeStrings;
		class PullOutputDataToMap : public OnDiskMatrixTransposer::PullOutputData {
		public:
			PullOutputDataToMap(const vector<string>* samples) : _samples(samples), _sampGts(new SampleToGenotypeStrings()) {}
			virtual ~PullOutputDataToMap() { delete _sampGts; }

			SampleToGenotypeStrings* getSampleGenotypes() const { return _sampGts; }

			virtual void setNextDataPoint(const string& data, int ind) {
				const string& samp = (*_samples)[ind];
				(*_sampGts)[samp] = data;
			}

		private:
			const vector<string>* _samples;
			SampleToGenotypeStrings* _sampGts;
		};
	};

}

#endif
