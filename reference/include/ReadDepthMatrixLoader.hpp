#ifndef __READ_DEPTH_MATRIX_LOADER_H__
#define __READ_DEPTH_MATRIX_LOADER_H__

#include "Interval.hpp"

#include <DataLoader.hpp>
#include <PreciseNonNegativeReal.hpp>
#include <UnivariateQuantitativeData.hpp>

#include <vector>
#include <set>
using namespace std;

namespace XHMM {

	class ReadDepthMatrixLoader : public HMM_PP::DataLoader<HMM_PP::UnivariateQuantitativeData> {

	public:
		ReadDepthMatrixLoader(string dataFile, const set<string>* keepIDs = NULL, double maxQuantitativeDataVal = HMM_PP::PreciseNonNegativeReal<double>::REAL_INFINITY);
		virtual ~ReadDepthMatrixLoader();

		virtual HMM_PP::UnivariateQuantitativeData* next();

		inline uint getNumTargets() const { return _targets->size(); }
		inline const Interval& getTarget(uint i) const { return (*_targets)[i]; }

		bool hasTargetsForChr(string chr) const;

		uint getChrStartTargetIndex(string chr) const;
		uint getChrStopTargetIndex(string chr) const;

		const set<uint>& getChrStopTargetIndices() const { return *_listAllChrStopInds; }

	protected:
		vector<Interval>* _targets;

		map<string, uint>* _chrStartInds;
		map<string, uint>* _chrStopInds;

		set<uint>* _listAllChrStopInds;

		virtual void setValFromStream(HMM_PP::UnivariateQuantitativeData::InputType& f, istream& stream);

	private:
		const double _maxQuantitativeDataVal;

		void readTargets();
		void printTargets(ostream& stream);
	};
}

#endif
