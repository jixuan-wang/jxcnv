#include "MultivariateQuantitativeData.hpp"

HMM_PP::MultivariateQuantitativeData::MultivariateQuantitativeData()
: Data() {
}

HMM_PP::MultivariateQuantitativeData::~MultivariateQuantitativeData() {
}

void HMM_PP::MultivariateQuantitativeData::clearObservations() {
	_mdata_qt.clear();
}

void HMM_PP::MultivariateQuantitativeData::setNumObservations(const uint nobs) {
	Data::setNumObservations(nobs);

	_mdata_qt.clear();
	_mdata_qt.resize(_n);
}

void HMM_PP::MultivariateQuantitativeData::setDatapoint(const uint i, const vector<double>& t) {
	if (i < _n) _mdata_qt[i] = t;
}

void HMM_PP::MultivariateQuantitativeData::addDatapoint(const vector<double>& t) {
	_mdata_qt.push_back(t);
	_n = _mdata_qt.size();
}

double HMM_PP::MultivariateQuantitativeData::qt(const uint t, const uint j) {
	return _mdata_qt[t][j];
}

void HMM_PP::MultivariateQuantitativeData::printDatapoint(ostream& stream, const uint t) const {
	const vector<double>& tData = _mdata_qt[t];
	for (vector<double>::const_iterator it = tData.begin(); it != tData.end(); ++it)
		stream << *it << '\t';
}
