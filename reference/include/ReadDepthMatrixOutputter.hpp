#ifndef __READ_DEPTH_MATRIX_OUTPUTTER_H__
#define __READ_DEPTH_MATRIX_OUTPUTTER_H__

#include <DataOutputter.hpp>

#include <string>
using namespace std;

namespace HMM_PP {
	class Data;
	class ostreamWriter;
}

namespace XHMM {
	class ReadDepthMatrixLoader;

	class ReadDepthMatrixOutputter : public HMM_PP::DataOutputter<HMM_PP::Data> {

	public:
		ReadDepthMatrixOutputter(string outFile, const ReadDepthMatrixLoader* readDepthLoader, int precision = 8);
		virtual ~ReadDepthMatrixOutputter();

		virtual void printDataType(const HMM_PP::Data* d);

		static void printTargets(ostream& stream, const ReadDepthMatrixLoader* readDepthLoader);

	protected:
		HMM_PP::ostreamWriter* _outStream;
		int _precision;
	};

}

#endif
