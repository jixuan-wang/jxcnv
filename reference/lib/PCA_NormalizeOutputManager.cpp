#include "PCA_NormalizeOutputManager.hpp"
#include "xhmmInputManager.hpp"
#include "xhmmOutputManager.hpp"
#include "ReadDepthMatrixLoader.hpp"
#include "ReadDepthMatrixOutputter.hpp"

#include <NamedMatrix.hpp>
#include <MatrixDecomp.hpp>
#include <PreciseNonNegativeReal.hpp>

#include <sstream>
#include <iomanip>
using namespace std;

const string XHMM::PCA_NormalizeOutputManager::PRINC_COMP_SUFFIX = ".PC.txt";
const string XHMM::PCA_NormalizeOutputManager::PC_SAMPLE_LOADINGS_SUFFIX = ".PC_LOADINGS.txt";
const string XHMM::PCA_NormalizeOutputManager::PC_SD_SUFFIX = ".PC_SD.txt";

XHMM::PCA_NormalizeOutputManager::PCA_NormalizeOutputManager(string PCAfile, string saveMemoryDir, int precision)
: _PCfile(PCAfile + PRINC_COMP_SUFFIX),
  _PCloadingsFile(PCAfile + PC_SAMPLE_LOADINGS_SUFFIX),
  _PCstdDevFile(PCAfile + PC_SD_SUFFIX),
  _saveMemoryDir(saveMemoryDir),
  _precision(precision) {
}

XHMM::PCA_NormalizeOutputManager::~PCA_NormalizeOutputManager() {
}

void XHMM::PCA_NormalizeOutputManager::PCA(string rdFile) const {
	xhmmInputManager::LoadedReadDepths loadRD = xhmmInputManager::loadReadDepthsFromFile(rdFile);
	delete loadRD._excludedTargets;
	delete loadRD._excludedSamples;
	HMM_PP::DoubleMat* rdMat = loadRD._rd;

	HMM_PP::ullintPair rdMatDims = rdMat->getDims();

	const uint numSamples = rdMatDims.first;
	vector<string>* samples = new vector<string>(numSamples);
	for (uint i = 0; i < numSamples; ++i)
		(*samples)[i] = rdMat->rowName(i);

	const uint numTargets = rdMatDims.second;
	vector<string>* targets = new vector<string>(numTargets);
	for (uint j = 0; j < numTargets; ++j)
		(*targets)[j] = rdMat->colName(j);

	if (numSamples > numTargets) {
		/* Otherwise, may get LAPACK error: "** On entry to DGESDD parameter number 12 had an illegal value"
		 *
		 * See:
		 * https://bugzilla.redhat.com/show_bug.cgi?id=169399
		 * where it claims that, for an m x n-dimensional matrix, this is an issue when m is sufficiently larger than n.
		 *
		 * TODO: Could in principle transpose the matrix here and run that way,
		 * but then the downstream use / semantic interpretation of the singular values, singular vectors, and loadings would need to be revised!!!
		 */
		stringstream str;
		str << "Do not want to run PCA+normalization (using SVD decomposition) for # samples [" << numSamples << "] > # targets [" << numTargets << "]";
		throw new Exception(str.str());
	}

	HMM_PP::MatrixDecomp::SVDecomp<double> decomp = HMM_PP::MatrixDecomp::svd<double>(rdMat, true, _saveMemoryDir); // get V-transpose as V
#if defined(DEBUG)
	HMM_PP::BaseRealMat* sigma = decomp.D->diag();
	HMM_PP::BaseRealMat* sigmaVtranspose = *sigma * *decomp.V;
	delete sigma;
	HMM_PP::BaseRealMat* testMat = *decomp.U * *sigmaVtranspose;
	delete sigmaVtranspose;
	HMM_PP::BaseRealMat* diff = *testMat - *rdMat;
	delete testMat;
	double maxDiff = max(abs(diff->min()), abs(diff->max()));
	delete diff;
	if (!HMM_PP::PreciseNonNegativeReal<double>::equalsReal(maxDiff, 0.0)) {
		stringstream str;
		str << "SVD decomposition does not reconstruct the original read-depth matrix: maxDiff= " << maxDiff;
		throw new Exception(str.str());
	}
#endif

	HMM_PP::DoubleMat* Dmat = decomp.D->asMatrix(true); // get matrix of column vector
	Dmat->setMatrixName("D");
	Dmat->setColName(0, "SD");

	// Set Principal Component (PC) titles for V, U, and D:
	for (uint pc = 0; pc < numSamples; ++pc) {
		stringstream PCstr;
		PCstr << "PC" << (pc + 1);
		decomp.V->setRowName(pc, PCstr.str());
		decomp.U->setColName(pc, PCstr.str());
		Dmat->setRowName(pc, PCstr.str());
	}

	// Format and output V-transpose:
	for (uint targ = 0; targ < targets->size(); ++targ)
		decomp.V->setColName(targ, (*targets)[targ]);
	delete targets;

	xhmmOutputManager::printMatrix(*(decomp.V), _PCfile, true, true, _precision);

	// Format and output D:
	xhmmOutputManager::printMatrix(*Dmat, _PCstdDevFile, true, true, _precision);
	delete Dmat;

	// Format and output U:
	for (uint i = 0; i < numSamples; ++i)
		decomp.U->setRowName(i, (*samples)[i]);
	delete samples;

	HMM_PP::DoubleMat* Utranspose = decomp.U->transpose();
	xhmmOutputManager::printMatrix(*Utranspose, _PCloadingsFile, true, true, _precision);
	delete Utranspose;

	decomp.deleteData();
}

void XHMM::PCA_NormalizeOutputManager::normalize(string rdFile, string outFile, const enum_PCnormalizeMethod& normMethod, const uint numPCtoRemove, const double PVE_mean_factor, const double PVE_contrib) const {
	uint removeFirstPC = 0;

	HMM_PP::DoubleMat* princCompStdDev = xhmmInputManager::MatrixReader<double>::readMatrixFromFile(_PCstdDevFile);
	const uint numPrincComp = princCompStdDev->nrow();

	if (princCompStdDev->ncol() != 1)
		throw new Exception("PC stdDev file " + _PCstdDevFile + " should contain a column vector");

	if (normMethod == PCnormalizeMethod_arg_numPCtoRemove) {
		if (numPCtoRemove <= 0)
			throw new Exception("Must specify a positive number of numPCtoRemove");
		removeFirstPC = numPCtoRemove;
	}
	else if (normMethod == PCnormalizeMethod_arg_PVE_mean || normMethod == PCnormalizeMethod_arg_PVE_contrib) {
		HMM_PP::NamedVector<real>* princCompVar = new HMM_PP::NamedVector<real>(numPrincComp);
		for (uint i = 0; i < numPrincComp; ++i)
			(*princCompVar)[i] = (*princCompStdDev)(i, 0);
		// Square the standard deviations to get the variances:
		(*princCompVar) *= (*princCompVar);

		real totalVariance = princCompVar->sum();

		if (normMethod == PCnormalizeMethod_arg_PVE_mean) {
			/* Concept follows (Everitt and Dunn, 2001) described at:
			 * http://public.lanl.gov/mewall/kluwer2002.html
			 */
			if (PVE_mean_factor <= 0)
				throw new Exception("Must specify positive PVE_mean_factor argument");

			real scaledMeanVariance = (totalVariance / numPrincComp) * PVE_mean_factor;
			stringstream scaledMeanVarStr;
			scaledMeanVarStr << setiosflags(ios::fixed) << setprecision(2) << scaledMeanVariance.getLog10Value();
			cerr << "Removing all components with variance >= " << PVE_mean_factor << " * mean(variance) = 10^" << scaledMeanVarStr.str() << endl;

			for (uint pc = 0; pc < numPrincComp; ++pc) {
				if ((*princCompVar)[pc] >= scaledMeanVariance)
					removeFirstPC = pc + 1; // # of PCs to remove: pc index + 1
				else
					// Since the principal components are sorted in decreasing order of standard deviations, once we go below the mean we stay there:
					break;
			}
		}
		else if (normMethod == PCnormalizeMethod_arg_PVE_contrib) {
			if (PVE_contrib <= 0 || PVE_contrib >= 100)
				throw new Exception("Must specify PVE_contrib argument in the open interval (0, 100)");
			real sumVarianceToAchieve = totalVariance * PVE_contrib / 100;

			real cumSumVariance = 0;
			for (uint pc = 0; pc < numPrincComp; ++pc) {
				cumSumVariance += (*princCompVar)[pc];
				if (cumSumVariance >= sumVarianceToAchieve) {
					removeFirstPC = pc + 1; // # of PCs to remove: pc index + 1
					break;
				}
			}

			stringstream PVEstr;
			PVEstr << setiosflags(ios::fixed) << setprecision(2) << (cumSumVariance * 100 / totalVariance);
			cerr << "Removing the smallest set of highest components that explain at least " << PVE_contrib << "% of the variance: " << PVEstr.str() << "%" << endl;
		}

		delete princCompVar;

		if (removeFirstPC == 0) {
			stringstream errStr;
			errStr << "The variance in the PCA matrix is 10^";
			errStr << setiosflags(ios::fixed) << setprecision(2) << totalVariance.getLog10Value();
			errStr << ", so cannot find any components to remove. This is probably due to using too few samples (N= " << numPrincComp << "). Please try again with at least ~10 samples...";
			throw new Exception(errStr.str());
		}
	}

	if (removeFirstPC == 0)
		throw new Exception("Must remove > 0 PC");
	removeFirstPC = min(removeFirstPC, numPrincComp);
	delete princCompStdDev;
	cerr << "Removing first " << removeFirstPC << " principal components from read-depth data" << endl;

	string numRemovedFile = outFile + ".num_removed_PC.txt";
	xhmmOutputManager::printBasicContainer(vector<uint>(1, removeFirstPC), numRemovedFile);

	const HMM_PP::DoubleMat* princCompAsRows = xhmmInputManager::MatrixReader<double>::readMatrixFromFile(_PCfile, true, true, removeFirstPC);
	if (princCompAsRows->nrow() != removeFirstPC) {
		stringstream str;
		str << "Unable to read " << removeFirstPC << " rows from " << _PCfile << ", read only " << princCompAsRows->nrow() << " rows";
		throw new Exception(str.str());
	}

	ReadDepthMatrixLoader* dataLoader = new ReadDepthMatrixLoader(rdFile);

	if (princCompAsRows->ncol() != dataLoader->getNumTargets()) {
		stringstream str;
		str << "Number of targets in " << rdFile << " does not match " << _PCfile << " - expected: " << princCompAsRows->ncol() << " found: " << dataLoader->getNumTargets();
		throw new Exception(str.str());
	}
	for (uint i = 0; i < dataLoader->getNumTargets(); ++i)
		if (dataLoader->getTarget(i) != Interval(princCompAsRows->colName(i)))
			throw new Exception("Target in " + rdFile + " does not match " + _PCfile + " - expected: " + princCompAsRows->colName(i) + " found: " + dataLoader->getTarget(i).intervalString());

	ReadDepthMatrixOutputter* outputter = new ReadDepthMatrixOutputter(outFile, dataLoader);

	while (dataLoader->hasNext()) {
		HMM_PP::UnivariateQuantitativeData* d = dataLoader->next();

		HMM_PP::DoubleVec& data = d->getNonConstValues();
		for (uint pc = 0; pc < removeFirstPC; ++pc) {
			const HMM_PP::DoubleVec& subtractComponent = (*princCompAsRows)[pc];

			// Calculate the dot-product of the component and the data:
			HMM_PP::DoubleVec* subtractComponentCopy = new HMM_PP::DoubleVec(subtractComponent);
			double dataLoading = (*subtractComponentCopy *= data).sum();
			delete subtractComponentCopy;

			// Scale the data by this loading:
			HMM_PP::DoubleVec* dataScoresInComp = new HMM_PP::DoubleVec(subtractComponent);
			*dataScoresInComp *= dataLoading;

			// Subtract, from the data, the projection of the data into this [orthogonal] component:
			data -= *dataScoresInComp;
			delete dataScoresInComp;
		}

		outputter->printDataType(d);
		delete d;
	}

	delete princCompAsRows;
	delete outputter;
	delete dataLoader;
}
