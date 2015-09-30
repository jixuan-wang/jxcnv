#ifndef __TARGET_RANGE_H__
#define __TARGET_RANGE_H__

#include <utils.hpp>

#include <string>
#include <sstream>
using namespace std;

namespace XHMM {

	class TargetRange {

	public:
		TargetRange(const uint t1, const uint t2) : _startTarg(t1), _stopTarg(t2) {}
		TargetRange(const HMM_PP::ullintPair& t1_t2) : _startTarg(t1_t2.first), _stopTarg(t1_t2.second) {}
		TargetRange(const TargetRange& other) : _startTarg(other._startTarg), _stopTarg(other._stopTarg) {}

		~TargetRange() {}

		uint getStartTargIndex() const { return _startTarg; }
		uint getStopTargIndex() const { return _stopTarg; }
		uint getNumTargets() const { return _stopTarg - _startTarg + 1; }

		string getOneBasedTargetIndexSpan() const {
			stringstream str;
			str << (_startTarg + 1) << ".." << (_stopTarg + 1);
			return str.str();
		}

	private:
		const uint _startTarg;
		const uint _stopTarg;
	};

}

#endif
