#include "UnivariateQuantitativeData.hpp"
#include "sqlwrap.hpp"

void HMM_PP::UnivariateQuantitativeData::bindDatapointToDB(SQL& db, sqlite3_stmt* stmt, const string& data, const uint t) {
	db.bind_double(stmt, data, val(t));
}

void HMM_PP::UnivariateQuantitativeData::addDatapointFromDB(SQL& db, sqlite3_stmt* stmt, const uint queryPos) {
	this->UnivariateData<double>::addDatapoint(db.get_double(stmt, queryPos));
}
