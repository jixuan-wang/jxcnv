#include "MultivariateCategoricalData.hpp"

HMM_PP::MultivariateCategoricalData::MultivariateCategoricalData()
: CategoricalData() {
}

HMM_PP::MultivariateCategoricalData::~MultivariateCategoricalData() {
}

void HMM_PP::MultivariateCategoricalData::clearObservations() {
	_mdata_int.clear();
}

void HMM_PP::MultivariateCategoricalData::setNumObservations(const uint nobs) {
	CategoricalData::setNumObservations(nobs);

	_mdata_int.clear();
	_mdata_int.resize(_n);
}

void HMM_PP::MultivariateCategoricalData::setDatapoint(const uint i, const vector<uint>& t) {
	if (i < _n) _mdata_int[i] = t;
}

void HMM_PP::MultivariateCategoricalData::setDatapoint(const uint i, const vector<string>& t) {
	if (i < _n) _mdata_int[i] = ensureCategoryToType(t);
}

void HMM_PP::MultivariateCategoricalData::addDatapoint(const vector<uint>& t) {
	_mdata_int.push_back(t);
	_n = _mdata_int.size();
}

void HMM_PP::MultivariateCategoricalData::addDatapoint(const vector<string>& t) {
	_mdata_int.push_back(ensureCategoryToType(t));
	_n = _mdata_int.size();
}

uint HMM_PP::MultivariateCategoricalData::obs(const uint t, const uint j) {
	return _mdata_int[t][j];
}

void HMM_PP::MultivariateCategoricalData::printDatapoint(ostream& stream, const uint t) const {
	const vector<uint>& tData = _mdata_int[t];
	for (vector<uint>::const_iterator it = tData.begin(); it != tData.end(); ++it)
		stream << CategoricalData::typeToCategory(*it) << '\t';
}
