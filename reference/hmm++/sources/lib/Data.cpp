#include "Data.hpp"
#include "Exception.hpp"

#include <sstream>
using namespace std;

HMM_PP::Data::Data()
: _n(0), _id("") {
}

HMM_PP::Data::~Data() {
}

void HMM_PP::Data::clear() {
	_n = 0;
	this->clearObservations();
}

uint HMM_PP::Data::checkStartStopRange(const uint t1, const vector<uint>& statePath) const {
	const uint t2 = t1 + statePath.size() - 1;
	if (t1 > _n-1 || t2 > _n-1) {
		stringstream str;
		str << "Invalid t1= " << t1 << ", t2= " << t2;
		throw new Exception(str.str());
	}

	return t2;
}

void HMM_PP::Data::checkStartStopRange(const uint t1, const uint t2) const {
	if (t1 > _n-1 || t2 > _n-1 || t2 < t1) {
		stringstream str;
		str << "Invalid t1= " << t1 << ", t2= " << t2;
		throw new Exception(str.str());
	}
}
