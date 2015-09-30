#include "UnivariateCategoricalData.hpp"
#include "sqlwrap.hpp"
#include "Exception.hpp"

#include <istream>
using namespace std;

HMM_PP::UnivariateCategoricalData::UnivariateCategoricalData()
: Data(), UnivariateData<uint>(), CategoricalData(), DeclareInputType<string>() {
}

HMM_PP::UnivariateCategoricalData::~UnivariateCategoricalData() {
}

void HMM_PP::UnivariateCategoricalData::setDatapoint(const uint i, const string& t) {
	this->UnivariateData<uint>::setDatapoint(i, categoryToType(t));
}

void HMM_PP::UnivariateCategoricalData::addDatapoint(const string& t) {
	this->UnivariateData<uint>::addDatapoint(categoryToType(t));
}

void HMM_PP::UnivariateCategoricalData::printDatapoint(ostream& stream, const uint t) const {
	stream << CategoricalData::typeToCategory(val(t));
}

void HMM_PP::UnivariateCategoricalData::bindDatapointToDB(SQL& db, sqlite3_stmt* stmt, const string& data, const uint t) {
	db.bind_int(stmt, data, val(t));
}

void HMM_PP::UnivariateCategoricalData::addDatapointFromDB(SQL& db, sqlite3_stmt* stmt, const uint queryPos) {
	this->UnivariateData<uint>::addDatapoint(db.get_int(stmt, queryPos));
}
