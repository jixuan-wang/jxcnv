#ifndef __DISTRIBUTION_STATISTICS_H__
#define __DISTRIBUTION_STATISTICS_H__

#include "params.hpp"

#include <cmath>
using namespace std;

namespace HMM_PP {

	template<class T>
	class DistributionStatistics {

	public:
		DistributionStatistics(bool calcStdDev = true);

		~DistributionStatistics() {}

		void observeVal(const T& val);

		T getSum() const { return _sumVals; }
		T getMean() const;
		T getStdDev(bool unbiased = true) const;

	private:
		uint _count;

		T _sumVals;

		bool _calcStdDev;
		T _sumOfSquareVals;
	};

	template<class T>
	DistributionStatistics<T>::DistributionStatistics(bool calcStdDev)
	: _count(0), _sumVals(0), _calcStdDev(calcStdDev), _sumOfSquareVals(0) {
	}

	template<class T>
	void DistributionStatistics<T>::observeVal(const T& val) {
		++_count;
		_sumVals += val;

		if (_calcStdDev)
			_sumOfSquareVals += (val * val);
	}

	template<class T>
	T DistributionStatistics<T>::getMean() const {
		T mean = 0.0;

		if (_count > 0)
			mean = _sumVals / _count;

		return mean;
	}

	template<class T>
	T DistributionStatistics<T>::getStdDev(bool unbiased) const {
		T stdDev = 0.0;

		if (_calcStdDev && _count > 0) {
			T normalizedSquareOfSum = (_sumVals * _sumVals) / _count;

			uint unBiasedCount = _count;
			if (unbiased && _count > 1)
				unBiasedCount = _count - 1;

			if (unBiasedCount > 0)
				stdDev = sqrt((_sumOfSquareVals - normalizedSquareOfSum) / unBiasedCount);
		}

		return stdDev;
	}
}

#endif
