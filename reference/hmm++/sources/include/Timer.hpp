#ifndef __TIMER_H__
#define __TIMER_H__

#include "params.hpp"

#include <ctime>
using namespace std;

#define CUR_TIME time(NULL)

namespace HMM_PP {

	class Timer {
	public:
		Timer() : _startTime(CUR_TIME) {}
		~Timer() {}

		uint getDuration() const {
			return CUR_TIME - _startTime;
		}

	private:
		time_t _startTime;
	};

}

#endif
