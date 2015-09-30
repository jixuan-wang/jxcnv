#include "Interval.hpp"

#include <params.hpp>
#include <PreciseNonNegativeReal.hpp>
#include <Exception.hpp>

#include <sstream>
#include <string>
#include <list>
#include <cmath>
#include <set>
#include <map>
using namespace std;

map<string, uint> XHMM::Interval::CHR_TO_INDEX = map<string, uint>();

set<string> XHMM::Interval::SPAN_DELIMS = createSpanDelimsSet();

const double XHMM::Interval::KB = 1000.0;

set<string> XHMM::Interval::createSpanDelimsSet() {
	set<string> spanDelims;
	spanDelims.insert(string(1, DEFAULT_SPAN_DELIM));
	spanDelims.insert("..");

	return spanDelims;
}

XHMM::Interval::Interval(string chr, const uint bp1, const uint bp2)
: _chr(chr), _bp1(bp1), _bp2(bp2) {
	checkSpan();
}

XHMM::Interval::Interval(const Interval& i)
: _chr(i._chr), _bp1(i._bp1), _bp2(i._bp2) {
}

XHMM::Interval::Interval(const string& chrBp1Bp2) {
	init(chrBp1Bp2);
}

void XHMM::Interval::init(const string& chrBp1Bp2) {
	string::size_type chrEndPlusOne = chrBp1Bp2.find(CHR_SPAN_DELIM);
	if (chrEndPlusOne == string::npos) {
		stringstream str;
		str << "Cannot find '" << CHR_SPAN_DELIM << "' in target string " << chrBp1Bp2;
		throw new Exception(str.str());
	}
	_chr = chrBp1Bp2.substr(0, chrEndPlusOne);
	string bpSpanStr = chrBp1Bp2.substr(chrEndPlusOne + 1); // 1 == length of CHR_SPAN_DELIM

	string::size_type bp1EndPlusOne = string::npos;
	string spanDelim;
	for (set<string>::const_iterator delimIt = SPAN_DELIMS.begin(); delimIt != SPAN_DELIMS.end(); ++delimIt) {
		const string& delim = *delimIt;
		bp1EndPlusOne = bpSpanStr.find(delim);
		if (bp1EndPlusOne != string::npos) {
			spanDelim = delim;
			break;
		}
	}

	if (bp1EndPlusOne == string::npos) {
		stringstream bpStream(bpSpanStr);
		bpStream >> _bp1;
		_bp2 = _bp1;
		if (!bpStream) {
			stringstream str;
			str << "Cannot parse single base pair from target string " << chrBp1Bp2;
			throw new Exception(str.str());
		}
	}
	else {
		stringstream bp1Stream(bpSpanStr.substr(0, bp1EndPlusOne));
		bp1Stream >> _bp1;
		if (!bp1Stream) {
			stringstream str;
			str << "Cannot parse start base pair from target string " << chrBp1Bp2;
			throw new Exception(str.str());
		}

		stringstream bp2Stream(bpSpanStr.substr(bp1EndPlusOne + spanDelim.size(), bpSpanStr.size()));
		bp2Stream >> _bp2;
		if (!bp2Stream) {
			stringstream str;
			str << "Cannot parse stop base pair from target string " << chrBp1Bp2;
			throw new Exception(str.str());
		}
	}

	checkSpan();
}

XHMM::Interval::Interval(stringstream& lineStream) {
	list<string>* lineTokens = new list<string>();

	while (lineStream && !lineStream.eof()) {
		string str;
		lineStream >> str;
		lineTokens->push_back(str);
	}

	if (lineTokens->empty())
		throw new Exception("Cannot parse Interval from empty line");
	else if (lineTokens->size() < 3) {
		string chrBp1Bp2 = *(lineTokens->begin());
		init(chrBp1Bp2);
	}
	else {
		list<string>::const_iterator lineIt = lineTokens->begin();
		string chr = *(lineIt++);
		string bp1Str = *(lineIt++);
		string bp2Str = *(lineIt++);

		init(chr + CHR_SPAN_DELIM + bp1Str + DEFAULT_SPAN_DELIM + bp2Str);
	}

	delete lineTokens;
}

XHMM::Interval::~Interval() {
}

string XHMM::Interval::intervalString() const {
	stringstream str;

	str << _chr << CHR_SPAN_DELIM << _bp1;
	if (_bp2 != _bp1)
		str << DEFAULT_SPAN_DELIM << _bp2;

	return str.str();
}

uint XHMM::Interval::span() const {
	return (_bp2 - _bp1 + 1);
}

BaseReal XHMM::Interval::spanKB() const {
	return span() / KB;
}

uint XHMM::Interval::midpoint() const {
	return (_bp1 + _bp2) / 2;
}

uint XHMM::Interval::distance(const Interval& other, bool abutIsZero) const {
	if (this->_chr != other._chr)
		return UINT_INFINITY;

	// Subtract 1, since "0" distance would be for directly abutting targets:
	if (this->_bp2 < other._bp1) { // this is before other
		uint dist = other._bp1 - this->_bp2;
		if (abutIsZero)
			--dist;
		return dist;
	}
	else if (other._bp2 < this->_bp1) { // other is before this
		uint dist = this->_bp1 - other._bp2;
		if (abutIsZero)
			--dist;
		return dist;
	}

	// The segments overlap (or abut):
	return 0;
}

XHMM::Interval XHMM::Interval::operator+(const Interval& other) const {
	if (this->_chr != other._chr) {
		stringstream str;
		str << "Cannot add on target " << other << " to " << *this;
		throw new Exception(str.str());
	}

	return Interval(this->_chr, min(this->_bp1, other._bp1), max(this->_bp2, other._bp2));
}

bool XHMM::Interval::operator==(const Interval& i) const {
	return (!(*this < i) && !(i < *this));
}

bool XHMM::Interval::operator!=(const Interval& i) const {
	return !(*this == i);
}

bool XHMM::Interval::operator<(const Interval& i) const {
	if (this->_chr != i._chr)
		return LESS_CHR_STRING(this->_chr, i._chr);

	if (this->_bp1 != i._bp1)
		return (this->_bp1 < i._bp1);

	return (this->_bp2 < i._bp2);
}

XHMM::Interval& XHMM::Interval::operator=(const Interval& i) {
	_chr = i._chr;
	_bp1 = i._bp1;
	_bp2 = i._bp2;

	return *this;
}

void XHMM::Interval::checkSpan() {
	if (_bp2 < _bp1) {
		stringstream str;
		str << "Target interval must start BEFORE it ends: " << *this;
		throw new Exception(str.str());
	}
}

bool XHMM::Interval::LESS_CHR_STRING(const string& chr1, const string& chr2) {
	map<string, uint>::const_iterator chr1It = CHR_TO_INDEX.find(chr1);
	map<string, uint>::const_iterator chr2It = CHR_TO_INDEX.find(chr2);

	if (chr1It != CHR_TO_INDEX.end()) {
		if (chr2It != CHR_TO_INDEX.end())
			return chr1It->second < chr2It->second;
		else
			return true; // indexed before non-indexed
	}
	else {
		if (chr2It != CHR_TO_INDEX.end())
			return false; // non-indexed after indexed
		else
			return LESS_CHR_NUMBER_STRING(chr1, chr2);
	}
}

bool XHMM::Interval::LESS_CHR_NUMBER_STRING(const string& chr1, const string& chr2) {
	stringstream s1;
	s1 << chr1;
	int chr1Number;
	s1 >> chr1Number;
	bool chr1IsNumber = !s1.fail();

	stringstream s2;
	s2 << chr2;
	int chr2Number;
	s2 >> chr2Number;
	bool chr2IsNumber = !s2.fail();

	if (chr1IsNumber) {
		if (chr2IsNumber)
			return chr1Number < chr2Number;
		else
			return true; // numbers come before non-numbers
	}
	else { // !chr1IsNumber
		if (chr2IsNumber)
			return false; // non-numbers come after numbers
		else
			return chr1 < chr2;
	}
}
