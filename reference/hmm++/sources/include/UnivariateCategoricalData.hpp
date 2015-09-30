#ifndef __UNIVARIATE_CATEGORICAL_DATA_H__
#define __UNIVARIATE_CATEGORICAL_DATA_H__

#include "params.hpp"
#include "UnivariateData.hpp"
#include "CategoricalData.hpp"
#include "DeclareInputType.hpp"

#include <vector>
#include <string>
using namespace std;

namespace HMM_PP {
	class Model;

	// Observations for a particular sequence
	class UnivariateCategoricalData : public UnivariateData<uint>, public CategoricalData, public DeclareInputType<string> {

	public:
		UnivariateCategoricalData();
		virtual ~UnivariateCategoricalData();

		void setDatapoint(const uint i, const string& t);
		void addDatapoint(const string& t);

		virtual void printDatapoint(ostream& stream, const uint t) const;

		virtual void bindDatapointToDB(SQL& db, sqlite3_stmt* stmt, const string& data, const uint t);
		virtual void addDatapointFromDB(SQL& db, sqlite3_stmt* stmt, const uint queryPos);
	};
}

#endif
