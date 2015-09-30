#include "CategoricalData.hpp"
#include "Exception.hpp"

map<string, uint> HMM_PP::CategoricalData::_catmap;
map<uint, string> HMM_PP::CategoricalData::_catrmap;

HMM_PP::CategoricalData::~CategoricalData() {
}

uint HMM_PP::CategoricalData::ensureCategoryToType(const string& c) {
	if (c == "")
		throw new Exception("Invalid empty category name");

	map<string, uint>::iterator i = _catmap.find(c);
	if (i == _catmap.end()) {
		uint tmp = _catmap.size();
		_catmap[c] = tmp;
		_catrmap[tmp] = c;
	}

	return _catmap[c];
}

uint HMM_PP::CategoricalData::categoryToType(const string& c) {
	map<string, uint>::const_iterator findIt = _catmap.find(c);
	if (findIt == _catmap.end())
		throw new Exception("Invalid categorical state : " + c);

	return findIt->second;
}

string HMM_PP::CategoricalData::typeToCategory(const uint i) {
	return _catrmap[i];
}

bool HMM_PP::CategoricalData::isValidCategory(const string& c) {
	return _catmap.find(c) != _catmap.end();
}

vector<uint> HMM_PP::CategoricalData::ensureCategoryToType(const vector<string>& c) {
	vector<uint> r(c.size());
	for (uint i = 0; i < r.size(); ++i) r[i] = ensureCategoryToType(c[i]);
	return r;
}

void HMM_PP::CategoricalData::clearCategories() {
	_catmap.clear();
	_catrmap.clear();
}
