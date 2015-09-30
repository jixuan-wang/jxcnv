#include "ReadDepthMatrixLoader.hpp"

#include <Exception.hpp>
#include <Data.hpp>
#include <utils.hpp>

#include <sstream>
#include <set>
using namespace std;

XHMM::ReadDepthMatrixLoader::ReadDepthMatrixLoader(string dataFile, const set<string>* keepIDs, double maxQuantitativeDataVal)
: HMM_PP::DataLoader<HMM_PP::UnivariateQuantitativeData>(dataFile, true, true, keepIDs), // read the targets header line, and the sample IDs
  _targets(new vector<Interval>()), _chrStartInds(new map<string, uint>()), _chrStopInds(new map<string, uint>()),
  _listAllChrStopInds(new set<uint>()),
  _maxQuantitativeDataVal(maxQuantitativeDataVal) {

	readTargets();
	if (_targets->empty())
		throw new Exception("Read 0 targets from matrix");
}

XHMM::ReadDepthMatrixLoader::~ReadDepthMatrixLoader() {
	delete _targets;
	delete _chrStartInds;
	delete _chrStopInds;
	delete _listAllChrStopInds;
}

void XHMM::ReadDepthMatrixLoader::readTargets() {
	stringstream* lineStream = new stringstream(*_header);
	delete _header;
	_header = NULL;

	Interval* prevTarget = NULL;
	while (*lineStream && !lineStream->eof()) {
		string targetString;
		*lineStream >> targetString;
		if (!*lineStream)
			throw new Exception("Data input stream failed while reading target " + targetString);

		Interval* curTarget = new Interval(targetString);
		if (_chrStopInds->find(curTarget->getChr()) != _chrStopInds->end())
			throw new Exception("MUST provide targets in order and GROUPED BY CHROMOSOME: target " + curTarget->intervalString() + " reverts to chromosome " + curTarget->getChr() + ", which was interrupted by targets in other chromosomes");
		_targets->push_back(*curTarget);

		if (prevTarget == NULL || curTarget->getChr() != prevTarget->getChr()) // curTarget is a new chromosome
			(*_chrStartInds)[curTarget->getChr()] = _targets->size() - 1;

		if (prevTarget != NULL) {
			if (curTarget->getChr() != prevTarget->getChr()) // switched to a new chromosome, so can never go back to prevTarget->getChr()
				(*_chrStopInds)[prevTarget->getChr()] = _targets->size() - 2; // subtract 2 (and not 1) since already added curTarget
			else if (!(curTarget->getBp1() > prevTarget->getBp2()))
				throw new Exception("MUST provide NON-OVERLAPPING targets IN ORDER: target " + curTarget->intervalString() + " either overlaps with or precedes previous target " + prevTarget->intervalString());

			delete prevTarget;
		}
		prevTarget = curTarget;
	}

	if (prevTarget != NULL) {
		(*_chrStopInds)[prevTarget->getChr()] = _targets->size() - 1; // since now prevTarget points to the LAST target, which ends its chromosome
		delete prevTarget;
	}

	for (map<string, uint>::const_iterator stopIt = _chrStopInds->begin(); stopIt != _chrStopInds->end(); ++stopIt)
		_listAllChrStopInds->insert(stopIt->second);

	delete lineStream;
}

void XHMM::ReadDepthMatrixLoader::printTargets(ostream& stream) {
	for (vector<Interval>::const_iterator it = _targets->begin(); it != _targets->end(); ++it)
		stream << "Target: " << it->getChr() << " " << it->getBp1() << " " << it->getBp2() << endl;
}

HMM_PP::UnivariateQuantitativeData* XHMM::ReadDepthMatrixLoader::next() {
	HMM_PP::UnivariateQuantitativeData* d = HMM_PP::DataLoader<HMM_PP::UnivariateQuantitativeData>::next();

	if (d->getNumObservations() != _targets->size()) {
		stringstream str;
		str << "Each line must contain " << _targets->size() << " read-depths, one for each exome target";
		throw new Exception(str.str());
	}
	return d;
}

void XHMM::ReadDepthMatrixLoader::setValFromStream(HMM_PP::UnivariateQuantitativeData::InputType& f, istream& stream) {
	HMM_PP::DataLoader<HMM_PP::UnivariateQuantitativeData>::setValFromStream(f, stream);
	HMM_PP::capAbsoluteValue(f, _maxQuantitativeDataVal);
}

bool XHMM::ReadDepthMatrixLoader::hasTargetsForChr(string chr) const {
	return _chrStartInds->find(chr) != _chrStartInds->end() && _chrStopInds->find(chr) != _chrStopInds->end();
}

uint XHMM::ReadDepthMatrixLoader::getChrStartTargetIndex(string chr) const {
	map<string, uint>::const_iterator it = _chrStartInds->find(chr);
	if (it == _chrStartInds->end())
		throw new Exception("No such targeted chromosome: " + chr);
	return it->second;
}

uint XHMM::ReadDepthMatrixLoader::getChrStopTargetIndex(string chr) const {
	map<string, uint>::const_iterator it = _chrStopInds->find(chr);
	if (it == _chrStopInds->end())
		throw new Exception("No such targeted chromosome: " + chr);
	return it->second;
}
