#include "PrepareTargetsManager.hpp"
#include "xhmmOutputManager.hpp"
#include "xhmmInputManager.hpp"

#include <Exception.hpp>
#include <utils.hpp>

#include <sstream>
#include <list>
#include <string>
using namespace std;

XHMM::PrepareTargetsManager::PrepareTargetsManager() {
}

XHMM::PrepareTargetsManager::~PrepareTargetsManager() {
}

void XHMM::PrepareTargetsManager::prepareTargets(list<string>* targetsFiles, string mergedTargetsOutput) {
	xhmmInputManager::IntervalSet* allTargets = new xhmmInputManager::IntervalSet();

	for (list<string>::const_iterator targFileIt = targetsFiles->begin(); targFileIt != targetsFiles->end(); ++targFileIt) {
		xhmmInputManager::IntervalSet* targsToAdd = xhmmInputManager::readIntervalsFromFile(*targFileIt);
		allTargets->insert(targsToAdd->begin(), targsToAdd->end());
		delete targsToAdd;
	}

	list<Interval>* finalTargets = new list<Interval>();
	Interval* waitingTarg = NULL;

	// Target intervals are grouped by chromosome and sorted by start position and then stop position (due to Interval::operator< used by IntervalSet):
	for (xhmmInputManager::IntervalSet::const_iterator targIt = allTargets->begin(); targIt != allTargets->end(); ++targIt) {
		const Interval& targ = *targIt;

		if (waitingTarg != NULL && waitingTarg->getChr() == targ.getChr() && waitingTarg->abutsOrOverlaps(targ)) {
			// Merge overlapping (or abutting) intervals:
			*waitingTarg = *waitingTarg + targ;
		}
		else { // Stop waiting (output the old interval and start a new one):
			if (waitingTarg != NULL) {
				finalTargets->push_back(*waitingTarg);
				delete waitingTarg;
				waitingTarg = NULL;
			}

			waitingTarg = new Interval(targ);
		}
	}

	if (waitingTarg != NULL) {
		finalTargets->push_back(*waitingTarg);
		delete waitingTarg;
		waitingTarg = NULL;
	}

	xhmmOutputManager::printBasicContainer(*finalTargets, mergedTargetsOutput);
	delete finalTargets;
}
