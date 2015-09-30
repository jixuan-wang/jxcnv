#include <Model.hpp>
#include <InferredDataFit.hpp>

#include "PosteriorOutputManager.hpp"
#include "ReadDepthMatrixLoader.hpp"
#include "AuxiliaryCNVoutputManager.hpp"
#include "ReadDepthMatrixOutputter.hpp"

#define POST_FILE_SUFFIX_1 ".posteriors."
#define POST_FILE_SUFFIX_2 ".txt"

XHMM::PosteriorOutputManager::PosteriorOutputManager(string posteriorFile, const vector<uint>& types, const HMM_PP::Model* model, const ReadDepthMatrixLoader* readDepthLoader)
: _CNtypes() {

	for (vector<uint>::const_iterator typeIt = types.begin(); typeIt != types.end(); ++typeIt) {
		const uint type = *typeIt;

		stringstream typeFile;
		typeFile << posteriorFile << POST_FILE_SUFFIX_1 << model->getModelParams()->state(type) << POST_FILE_SUFFIX_2;
		HMM_PP::ostreamWriter* typeStream = HMM_PP::utils::getOstreamWriterFromFile(typeFile.str());

		ostream& stream = (*typeStream)();

		stream << setiosflags(ios::fixed);
		stream << setprecision(AuxiliaryCNVoutputManager::POSTERIOR_AND_RD_PRECISION);

		ReadDepthMatrixOutputter::printTargets(stream, readDepthLoader);

		_CNtypes.push_back(CNtypeAndWriter(type, typeStream));
	}
}

XHMM::PosteriorOutputManager::~PosteriorOutputManager() {
	for (list<CNtypeAndWriter>::iterator typeIt = _CNtypes.begin(); typeIt != _CNtypes.end(); ++typeIt)
		delete typeIt->second;
}

void XHMM::PosteriorOutputManager::printSamplePosteriors(const HMM_PP::InferredDataFit* idf) {
	for (list<CNtypeAndWriter>::iterator typeIt = _CNtypes.begin(); typeIt != _CNtypes.end(); ++typeIt) {
		const uint type = typeIt->first;
		HMM_PP::ostreamWriter* typeStream = typeIt->second;

		ostream& outStream = (*typeStream)();

		outStream
		<< idf->getData()->getId();

		const uint n = idf->getData()->getNumObservations();
		for (uint i = 0; i < n; ++i)
			outStream
			<< '\t'
			<< idf->getGamma(i, type);

		outStream << endl;
	}
}
