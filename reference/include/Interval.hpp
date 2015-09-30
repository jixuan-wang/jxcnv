#ifndef __INTERVAL_H__
#define __INTERVAL_H__

#include <params.hpp>

#include <string>
#include <map>
#include <sstream>
#include <set>
using namespace std;

namespace XHMM {

	class Interval {

	public:
		Interval(string chr, const uint bp1, const uint bp2);
		Interval(const Interval& i);
		Interval(const string& chrBp1Bp2);
		Interval(stringstream& lineStream);

		~Interval();

		static const double KB;

		const string& getChr() const { return _chr; }
		uint getBp1() const { return _bp1; }
		uint getBp2() const { return _bp2; }

		string intervalString() const;
		uint span() const;
		BaseReal spanKB() const;
		uint midpoint() const;

		/*
		 * If abutIsZero, then Intervals [1-2] and [3-4] will have
		 * a distance of 0 and so on for larger distances; otherwise, their distance is 1 and so on.
		 */
		uint distance(const Interval& other, bool abutIsZero = true) const;

		inline bool abutsOrOverlaps(const Interval& other) const { return this->distance(other) == 0; }
		inline bool strictlyOverlaps(const Interval& other) const { return this->distance(other, false) == 0; }

		// Return an Interval resulting from "adding" other onto this [that is, the smallest interval containing both this and other]:
		Interval operator+(const Interval& other) const;

		bool operator==(const Interval& i) const;
		bool operator!=(const Interval& i) const;

		bool operator<(const Interval& i) const;

		Interval& operator=(const Interval& i);

		static map<string, uint> CHR_TO_INDEX;

		inline friend ostream& operator<<(ostream& stream, const Interval& i) {
			return stream << i.intervalString();
		}

	private:
		string _chr;
		uint _bp1;
		uint _bp2;

		void init(const string& chrBp1Bp2);

		void checkSpan();

		static const char CHR_SPAN_DELIM = ':';

		static const char DEFAULT_SPAN_DELIM = '-';
		static set<string> SPAN_DELIMS;
		static set<string> createSpanDelimsSet();

		static bool LESS_CHR_STRING(const string& chr1, const string& chr2);
		static bool LESS_CHR_NUMBER_STRING(const string& chr1, const string& chr2);
	};
}

#endif
