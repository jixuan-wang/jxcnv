#ifndef __XHMM_OUTPUT_MANAGER_H__
#define __XHMM_OUTPUT_MANAGER_H__

#include <params.hpp>
#include <NamedMatrix.hpp>
#include <NamedVector.hpp>

#include <ostream>
#include <string>
using namespace std;

namespace XHMM {
	class CNVmodelParams;

	class xhmmOutputManager {

	public:
		static void printTransitionMatricesForDistancesFile(const CNVmodelParams* cnvParams, string distFile, ostream& out, int precision = 4);
		static void printInitProbs(const CNVmodelParams* cnvParams, ostream& out, int precision = 4);

		template<class RealType>
		static void printMatrix(const HMM_PP::NamedMatrix<RealType>& mat, string outFile, bool printRowNames = true, bool printColNames = true, int precision = 8);

		template<class RealType>
		static void printVector(const HMM_PP::NamedVector<RealType>& vec, string outFile, bool printNames = true, bool asColumn = false, int precision = 8);

		static void printWarning(const string& warn);

		template<class Container>
		static void printBasicContainer(const Container& cont, string outFile);
	};

	template<class RealType>
	void XHMM::xhmmOutputManager::printMatrix(const HMM_PP::NamedMatrix<RealType>& mat, string outFile, bool printRowNames, bool printColNames, int precision) {
		HMM_PP::ostreamWriter* outStream = HMM_PP::utils::getOstreamWriterFromFile(outFile);
		if (outStream == NULL)
			return;

		ostream& stream = (*outStream)();

		stream << setiosflags(ios::fixed);
		stream << setprecision(precision);

		mat.printDelimited(stream, printRowNames, printColNames);

		delete outStream;
	}

	template<class RealType>
	void XHMM::xhmmOutputManager::printVector(const HMM_PP::NamedVector<RealType>& vec, string outFile, bool printNames, bool asColumn, int precision) {
		HMM_PP::ostreamWriter* outStream = HMM_PP::utils::getOstreamWriterFromFile(outFile);
		if (outStream == NULL)
			return;

		ostream& stream = (*outStream)();

		stream << setiosflags(ios::fixed);
		stream << setprecision(precision);

		char delim = '\t';
		if (asColumn)
			delim = '\n';
		vec.printDelimited(stream, printNames, delim) << endl;

		delete outStream;
	}

	template<class Container>
	void XHMM::xhmmOutputManager::printBasicContainer(const Container& cont, string outFile) {
		HMM_PP::ostreamWriter* outStream = HMM_PP::utils::getOstreamWriterFromFile(outFile);
		if (outStream == NULL)
			return;

		ostream& stream = (*outStream)();
		for (typename Container::const_iterator it = cont.begin(); it != cont.end(); ++it)
			stream << *it << endl;

		delete outStream;
	}
}

#endif
