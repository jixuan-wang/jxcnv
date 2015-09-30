#ifndef __AUXILIARY_CNV_OUTPUT_MANAGER_H__
#define __AUXILIARY_CNV_OUTPUT_MANAGER_H__

#include "ReadDepthMatrixLoader.hpp"

#include <params.hpp>
#include <PreciseNonNegativeReal.hpp>
#include <utils.hpp>

namespace HMM_PP {
	class Model;
}

namespace XHMM {
	class CNVmodelParams;
	class SampleGenotypeQualsCalculator;

	class AuxiliaryCNVoutputManager {

	public:
		static const int POSTERIOR_AND_RD_PRECISION = 2;

		AuxiliaryCNVoutputManager
		(string auxOutputFile,
				uint auxUpstreamPrintTargs = 0, uint auxDownstreamPrintTargs = 0,
				bool printOrigData = false);

		~AuxiliaryCNVoutputManager();

		void setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader);

		void printCNVtargets(SampleGenotypeQualsCalculator* genotyper, const uint t1, const uint t2, const uint type, HMM_PP::Data* origData);
		void printCNVtargets(SampleGenotypeQualsCalculator* genotyper, const uint t1, const uint t2, const vector<uint>& types, HMM_PP::Data* origData);

	private:
		HMM_PP::ostreamWriter* _auxStream;

		uint _auxUpstreamPrintTargs;
		uint _auxDownstreamPrintTargs;

		HMM_PP::Model* _model;
		const ReadDepthMatrixLoader* _dataLoader;

		bool _printOrigData;

		void printTargetHeader();
	};

}

#endif
