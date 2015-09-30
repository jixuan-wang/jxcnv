#ifndef __COPY_NUMBER_VARIANT_H__
#define __COPY_NUMBER_VARIANT_H__

#include "Interval.hpp"
#include "TargetRange.hpp"
#include "Genotype.hpp"

#include <params.hpp>
#include <utils.hpp>

#include <map>
#include <vector>
#include <string>
using namespace std;

namespace XHMM {
	class ReadDepthMatrixLoader;

	class CopyNumberVariant {

	public:
		CopyNumberVariant(const vector<uint>& nonRefTypes, const HMM_PP::ullintPair* t1_t2, const ReadDepthMatrixLoader* dataLoader, Genotype::RealThresh genotypingThreshold, bool storeGenotypes = true, set<string>* discoveredSamples = NULL);
		~CopyNumberVariant();

		const Interval* getMergedTargets() const { return _mergedTargets; }
		const TargetRange* getTargetRange() const { return _targRange; }

		const Interval* getPreviousTarget() const { return _prevTarget; }
		const Interval* getPostTarget() const { return _postTarget; }

		// NOTE: No error-checking as to the Genotype added:
		void addGenotype(string sample, Genotype* gt);
		const Genotype* getGenotype(string sample) const;

		uint getNonRefTypeIndex(const uint type) const;

		uint getNonRefTypeCount(const uint type) const;
		uint getCalledSampleCount() const { return _calledSampleCount; }

		const Genotype::RealThresh& getGenotypingThreshold() const { return _genotypingThreshold; }

		const set<string>* getDiscoveredSamples() const { return _discoveredSamples; }

	private:
		Interval* _mergedTargets;
		Interval* _prevTarget;
		Interval* _postTarget;
		TargetRange* _targRange;

		map<uint, uint>* _nonRefTypeToOrderIndex;

		map<string, Genotype*>* _sampleGenotypes;

		map<uint, uint>* _nonRefTypeToSampleCount;
		uint _calledSampleCount;

		Genotype::RealThresh _genotypingThreshold;

		set<string>* _discoveredSamples;
	};

}

#endif
