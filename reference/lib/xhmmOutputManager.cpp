#include "xhmmOutputManager.hpp"
#include "Interval.hpp"

#include <CNVmodelParams.hpp>

void XHMM::xhmmOutputManager::printTransitionMatricesForDistancesFile(const CNVmodelParams* cnvParams, string distFile, ostream& out, int precision) {
	HMM_PP::istreamLineReader* in = HMM_PP::utils::getIstreamLineReaderFromFile(distFile);
	if (in == NULL)
		throw new Exception("Unable to read user-input distances");

	out << setprecision(precision);

	while (true) {
		out << "Enter inter-target distance (in KB) to calculate HMM transition matrix (CTRL-D to exit):" << endl;
		if (in->eof())
			break;

		string line;
		*in >> line;
		stringstream lineStream(line);

		double dist;
		lineStream >> dist;
		if (!lineStream)
			throw new Exception("Unable to read distance from input stream");
		if (dist < 0)
			throw new Exception("Cannot enter negative distances");

		HMM_PP::NamedMatrix<real>* mat = cnvParams->transitionMatrix(static_cast<uint>(dist * Interval::KB));

		out
		<< "HMM transition matrix for targets separated by " << dist << " KB (rows are target 't', columns are target 't+1'):" << endl
		<< *mat << endl;

		delete mat;
	}

	delete in;
}

void XHMM::xhmmOutputManager::printInitProbs(const CNVmodelParams* cnvParams, ostream& out, int precision) {
	HMM_PP::BaseRealVec printInitProbs = cnvParams->getInitProbs();
	cnvParams->giveStateNames(printInitProbs);

	out
	<< setprecision(precision)
	<< "Initial HMM probabilities:" << endl
	<< printInitProbs << endl << endl;
}

void XHMM::xhmmOutputManager::printWarning(const string& warn) {
	cerr << "WARNING: " << warn << endl;
}
