#ifndef __MULTIVARIATE_CATEGORICAL_DATA_H__
#define __MULTIVARIATE_CATEGORICAL_DATA_H__

#include "params.hpp"
#include "CategoricalData.hpp"

#include <vector>
using namespace std;

namespace HMM_PP {
	class Model;

	// Observations for a particular sequence
	class MultivariateCategoricalData : public CategoricalData {

	public:
		MultivariateCategoricalData();
		virtual ~MultivariateCategoricalData();

	protected:
		virtual void clearObservations();

	public:
		virtual void setNumObservations(const uint nobs);

		void setDatapoint(const uint i, const vector<uint>& t);
		void setDatapoint(const uint i, const vector<string>& t);

		void addDatapoint(const vector<uint>& t);
		void addDatapoint(const vector<string>& t);

		// get datapoint at time 't'
		uint obs(const uint t, const uint j);

		virtual void printDatapoint(ostream& stream, const uint t) const;

	private:
		vector<vector<uint> > _mdata_int;
	};
}

#endif
