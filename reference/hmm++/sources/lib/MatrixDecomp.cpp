#include "MatrixDecomp.hpp"
#include "Exception.hpp"
#include "VectorOnDisk.hpp"

#include <signal.h>

#include <algorithm>
#include <string>
#include <cctype>
using namespace std;

set<HMM_PP::LAPACKvector<double>*> HMM_PP::MatrixDecomp::_mappedMem;
void HMM_PP::MatrixDecomp::initSignalHandlers() {
	signal(SIGTERM, HMM_PP::MatrixDecomp::freeMappedMemory);
	//signal(SIGKILL, HMM_PP::MatrixDecomp::freeMappedMemory);
	signal(SIGSEGV, HMM_PP::MatrixDecomp::freeMappedMemory);
	signal(SIGCHLD, HMM_PP::MatrixDecomp::freeMappedMemory);
	signal(SIGPIPE, HMM_PP::MatrixDecomp::freeMappedMemory);
	signal(SIGINT,  HMM_PP::MatrixDecomp::freeMappedMemory);
	signal(SIGABRT, HMM_PP::MatrixDecomp::freeMappedMemory);
	signal(SIGALRM, HMM_PP::MatrixDecomp::freeMappedMemory);
}

void HMM_PP::MatrixDecomp::freeMappedMemory(int sig) {
	bool freedMappedMemory = false;
	for (set<LAPACKvector<double>*>::iterator i = _mappedMem.begin(); i != _mappedMem.end(); ++i) {
		freedMappedMemory = true;
		LAPACKvector<double>* mappedVec = *i;
		_mappedMem.erase(i);
		delete mappedVec;
	}
	_mappedMem.clear();

	if (!freedMappedMemory) {
		signal(sig, SIG_DFL);
		kill(getpid(), sig);
		return;
	}

	// Reset the signal handler:
	signal(sig, HMM_PP::MatrixDecomp::freeMappedMemory);

	stringstream str;
	str << "Exiting due to caught ";

	char* sigCstr = strsignal(sig);
	if (sigCstr) {
		string sigStr = sigCstr;
		transform(sigStr.begin(), sigStr.end(), sigStr.begin(), (int (*)(int)) std::toupper);
		str << sigStr;
	}
	else
		str << sig;

	str << " signal after memory mapping";
	throw new Exception(str.str());
}

HMM_PP::LAPACKvector<double>* HMM_PP::MatrixDecomp::allocateVector(ullint size, const string workDir) {
	HMM_PP::LAPACKvector<double>* retVec = NULL;

	if (workDir == "")
		retVec = new HMM_PP::StandardVector<double>(size);
	else {
		retVec = new HMM_PP::VectorOnDisk<double>(size, workDir);
		_mappedMem.insert(retVec);
	}

	return retVec;
}

void HMM_PP::MatrixDecomp::deleteAllocatedVector(HMM_PP::LAPACKvector<double>* vec) {
	_mappedMem.erase(vec);
	delete vec;
}
