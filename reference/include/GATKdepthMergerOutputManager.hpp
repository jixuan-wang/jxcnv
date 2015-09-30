#ifndef __GATK_DEPTH_MERGER_OUTPUT_MANAGER_H__
#define __GATK_DEPTH_MERGER_OUTPUT_MANAGER_H__

#include "Interval.hpp"
#include "xhmmInputManager.hpp"
#include "AuxiliaryCNVoutputManager.hpp"
#include "OnDiskMatrixTransposer.hpp"

#include <params.hpp>

#include <string>
#include <list>
#include <utility>
#include <vector>
using namespace std;

namespace HMM_PP {
	class istreamLineReader;
}

namespace XHMM {

	class GATKdepthMergerOutputManager {

	public:
		GATKdepthMergerOutputManager(list<string>* GATKdepthFiles, list<string>* GATKdepthListFiles, xhmmInputManager::OneToOneStrMap* mapSamples = NULL, string columnSuffix = "_mean_cvg");
		~GATKdepthMergerOutputManager();

		void mergeGATKdepths(string outRdFile, bool outputSamplesByTargets = true, int rdPrecision = AuxiliaryCNVoutputManager::POSTERIOR_AND_RD_PRECISION, bool transposeInMemory = true);

	private:
		const string _columnSuffix;

		list<Interval>* _allTargets;
		vector<string>* _samples;

		class GATKdepthReader;
		list<GATKdepthReader*>* _readers;

		class PrintOutputData : public OnDiskMatrixTransposer::PullOutputData {
		public:
			PrintOutputData(ostream& stream) : _stream(stream) {}
			virtual ~PrintOutputData() {}

			virtual void setNextDataPoint(const string& data, int ind) {
				_stream << '\t' << data;
			}

		private:
			ostream& _stream;
		};

		class GATKdepthReader {
		public:
			GATKdepthReader(string GATKdepthFile, string columnSuffix, bool verbose = true);
			~GATKdepthReader();

			typedef pair<string, double> SampleValue;
			typedef list<SampleValue> SampleValues;
			typedef pair<Interval, SampleValues*> TargetSampleValues;

			bool hasNextTarget();
			TargetSampleValues nextTargetValues(bool extractValues = true);

			list<Interval>* getRemainingTargets();

			typedef pair<string, uint> SampleIndex;
			const list<SampleIndex>* getSamples() const { return _samples; }

		private:
			HMM_PP::istreamLineReader* _stream;
			list<SampleIndex>* _samples;

			const string _columnSuffix;
			const uint _columnSuffixLength;

			void readSamplesFromHeader();
		};

	};

}

#endif
