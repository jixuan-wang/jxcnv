#ifndef __MULTIVARIATE_QUANTITATIVE_DATA_H__
#define __MULTIVARIATE_QUANTITATIVE_DATA_H__

#include "params.hpp"
#include "Data.hpp"

#include <vector>
using namespace std;

namespace HMM_PP {
	class Model;

	// Observations for a particular sequence
	class MultivariateQuantitativeData : public Data {

	public:
		MultivariateQuantitativeData();
		virtual ~MultivariateQuantitativeData();

	protected:
		virtual void clearObservations();

	public:
		virtual void setNumObservations(const uint nobs);

		void setDatapoint(const uint i, const vector<double>& t);

		void addDatapoint(const vector<double>& t);

		// get datapoint at time 't'
		double qt(const uint t, const uint j);

		virtual void printDatapoint(ostream& stream, const uint t) const;

	private:
		vector<vector<double> > _mdata_qt;
	};
}

#endif
