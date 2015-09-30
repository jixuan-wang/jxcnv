#ifndef __PREPARE_TARGETS_MANAGER_H__
#define __PREPARE_TARGETS_MANAGER_H__

#include "Interval.hpp"
#include "xhmmInputManager.hpp"

#include <params.hpp>

namespace XHMM {

	class PrepareTargetsManager {

	public:
		PrepareTargetsManager();
		~PrepareTargetsManager();

		void prepareTargets(list<string>* targetsFiles, string mergedTargetsOutput);

	private:
	};

}

#endif
