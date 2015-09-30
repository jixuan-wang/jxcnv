#ifndef __UNIVARIATE_DATA_NO_DB_H__
#define __UNIVARIATE_DATA_NO_DB_H__

#include "UnivariateData.hpp"
#include "DeclareInputType.hpp"

namespace HMM_PP {

	// Observations for a particular sequence
	template<class BaseType>
	class UnivariateDataNoDB : public UnivariateData<BaseType>, public DeclareInputType<BaseType> {

	public:
		UnivariateDataNoDB() : Data(), UnivariateData<BaseType>(), DeclareInputType<BaseType>() {}
		virtual ~UnivariateDataNoDB() {}

		virtual void bindDatapointToDB(SQL& db, sqlite3_stmt* stmt, const string& data, const uint t) { throw new Exception("DB functions not implemented"); }
		virtual void addDatapointFromDB(SQL& db, sqlite3_stmt* stmt, const uint queryPos) { throw new Exception("DB functions not implemented"); }
	};
}

#endif
