#ifndef __UNIVARIATE_QUANTITATIVE_DATA_H__
#define __UNIVARIATE_QUANTITATIVE_DATA_H__

#include "UnivariateDataNoDB.hpp"

namespace HMM_PP {

	// Observations for a particular sequence
	class UnivariateQuantitativeData : public UnivariateDataNoDB<double> {

	public:
		UnivariateQuantitativeData() : Data(), UnivariateDataNoDB<double>() {}
		virtual ~UnivariateQuantitativeData() {}

		virtual void bindDatapointToDB(SQL& db, sqlite3_stmt* stmt, const string& data, const uint t);
		virtual void addDatapointFromDB(SQL& db, sqlite3_stmt* stmt, const uint queryPos);
	};
}

#endif
