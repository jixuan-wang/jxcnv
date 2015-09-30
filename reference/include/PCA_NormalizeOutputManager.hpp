#ifndef __PCA_NORMALIZE_OUTPUT_MANAGER_H__
#define __PCA_NORMALIZE_OUTPUT_MANAGER_H__

#include "xhmmCmdline.h"

#include <params.hpp>

#include <string>
using namespace std;

namespace XHMM {

	class PCA_NormalizeOutputManager {

	public:
		PCA_NormalizeOutputManager(string PCAfile, string saveMemoryDir = "", int precision = 6);
		~PCA_NormalizeOutputManager();

		void PCA(string rdFile) const;
		void normalize(string rdFile, string outFile, const enum_PCnormalizeMethod& normMethod, const uint numPCtoRemove, const double PVE_mean_factor, const double PVE_contrib) const;

	private:
		const string _PCfile;
		const string _PCloadingsFile;
		const string _PCstdDevFile;

		const string _saveMemoryDir;

		const int _precision;

		static const string PRINC_COMP_SUFFIX;
		static const string PC_SAMPLE_LOADINGS_SUFFIX;
		static const string PC_SD_SUFFIX;
	};

}

#endif
