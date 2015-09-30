#include "ReadDepthMatrixOutputter.hpp"
#include "ReadDepthMatrixLoader.hpp"
#include "DiscoverOutputManager.hpp"

#include <Data.hpp>
#include <utils.hpp>

#include <iomanip>
using namespace std;

XHMM::ReadDepthMatrixOutputter::ReadDepthMatrixOutputter(string outFile, const ReadDepthMatrixLoader* readDepthLoader, int precision)
: _precision(precision) {

	_outStream = HMM_PP::utils::getOstreamWriterFromFile(outFile);
	ostream& stream = (*_outStream)();

	stream << setiosflags(ios::fixed);
	stream << setprecision(_precision);

	printTargets(stream, readDepthLoader);
}

XHMM::ReadDepthMatrixOutputter::~ReadDepthMatrixOutputter() {
	delete _outStream;
}

void XHMM::ReadDepthMatrixOutputter::printTargets(ostream& stream, const ReadDepthMatrixLoader* readDepthLoader) {
	stream << XHMM::DiscoverOutputManager::INTERVAL;

	for (uint i = 0; i < readDepthLoader->getNumTargets(); ++i)
		stream << '\t' << readDepthLoader->getTarget(i);

	stream << endl;
}

void XHMM::ReadDepthMatrixOutputter::printDataType(const HMM_PP::Data* d) {
	ostream& stream = (*_outStream)();

	stream << d->getId();

	const uint n = d->getNumObservations();
	for (uint i = 0; i < n; ++i) {
		stream << '\t';
		d->printDatapoint(stream, i);
	}

	stream << endl;
}
