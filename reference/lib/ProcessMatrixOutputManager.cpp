#include "ProcessMatrixOutputManager.hpp"
#include "xhmmInputManager.hpp"
#include "xhmmOutputManager.hpp"

#include <NamedMatrix.hpp>
#include <DistributionStatistics.hpp>
#include <utils.hpp>

XHMM::ProcessMatrixOutputManager::ProcessMatrixOutputManager() {
}

XHMM::ProcessMatrixOutputManager::~ProcessMatrixOutputManager() {
}

void XHMM::ProcessMatrixOutputManager::processReadDepthMatrix(string rdFile, string outFile,
		list<string>* excludeTargetsFiles, list<string>* excludeChromosomeTargets, list<string>* excludeSamplesFiles,
		const uint minTargetSize, const uint maxTargetSize,
		const double minMeanTargetRD, const double maxMeanTargetRD,
		const double minSdTargetRD, const double maxSdTargetRD,
		const double minMeanSampleRD, const double maxMeanSampleRD,
		const double minSdSampleRD, const double maxSdSampleRD,
		bool scaleDataBySum, const enum_scaleDataBySumType& scaleDataBySumType, const double scaleDataBySumFactor,
		bool log10Data, const double pseudoCountValForLog10Input,
		bool centerData, const enum_centerType& centerType,
		bool zScoreData,
		string outputExcludedTargets, string outputExcludedSamples) {

	// Combine all targets to exclude:
	xhmmInputManager::IntervalSet* excludeTargets = new xhmmInputManager::IntervalSet();
	for (list<string>::const_iterator targIt = excludeTargetsFiles->begin(); targIt != excludeTargetsFiles->end(); ++targIt) {
		xhmmInputManager::IntervalSet* exclude = xhmmInputManager::readIntervalsFromFile(*targIt);
		excludeTargets->insert(exclude->begin(), exclude->end());
		delete exclude;
	}
	delete excludeTargetsFiles;

	// Combine all target chromosomes to exclude:
	xhmmInputManager::StringSet* excludeTargetChromosomes = new xhmmInputManager::StringSet(excludeChromosomeTargets->begin(), excludeChromosomeTargets->end());
	delete excludeChromosomeTargets;

	// Combine all samples to exclude:
	xhmmInputManager::StringSet* excludeSamples = new xhmmInputManager::StringSet();
	for (list<string>::const_iterator sampIt = excludeSamplesFiles->begin(); sampIt != excludeSamplesFiles->end(); ++sampIt) {
		list<string>* exclude = xhmmInputManager::readStringsFromFile(*sampIt);
		excludeSamples->insert(exclude->begin(), exclude->end());
		delete exclude;
	}
	delete excludeSamplesFiles;

	// Read the initial matrix:
	xhmmInputManager::LoadedReadDepths loadRD = xhmmInputManager::loadReadDepthsFromFile(rdFile, excludeTargets, excludeTargetChromosomes, excludeSamples, minTargetSize, maxTargetSize);
	HMM_PP::DoubleMat* rdMat = loadRD._rd;
	xhmmInputManager::IntervalSet* excludedTargets = loadRD._excludedTargets;
	xhmmInputManager::StringSet* excludedSamples = loadRD._excludedSamples;

	delete excludeTargets;
	delete excludeTargetChromosomes;
	delete excludeSamples;

	rdMat = filterTargetProperties(rdMat,
			minMeanTargetRD, maxMeanTargetRD,
			minSdTargetRD, maxSdTargetRD,
			excludedTargets);

	rdMat = filterSampleProperties(rdMat,
			minMeanSampleRD, maxMeanSampleRD,
			minSdSampleRD, maxSdSampleRD,
			excludedSamples);

	if (scaleDataBySum) {
		switch (scaleDataBySumType) {
			case scaleDataBySumType_arg_target :{
				*rdMat *= scaleDataBySumFactor;
				rdMat->scaleByColumnSums();
				cerr << "Scaled by target sums * " << scaleDataBySumFactor << "." << endl;
				break;
			}
			case scaleDataBySumType_arg_sample :{
				*rdMat *= scaleDataBySumFactor;
				rdMat->scaleByRowSums();
				cerr << "Scaled by sample sums * " << scaleDataBySumFactor << "." << endl;
				break;
			}
			default :{
				XHMM::xhmmOutputManager::printWarning("No scaling performed since 'scaleDataBySumType' not given");
				break;
			}
		}
	}

	if (log10Data) {
		rdMat->log10(pseudoCountValForLog10Input);
		cerr << "Applied log10(max(values, 0) + " << pseudoCountValForLog10Input << ") to matrix values." << endl;
	}

	if (centerData) {
		switch (centerType) {
			case centerType_arg_target :{
				rdMat->centerColumnsByMean(zScoreData);
				cerr << "Centered by target mean" << (zScoreData ? " and calculated z-scores" : "") << "." << endl;
				break;
			}
			case centerType_arg_sample :{
				rdMat->centerRowsByMean(zScoreData);
				cerr << "Centered by sample mean" << (zScoreData ? " and calculated z-scores" : "") << "." << endl;
				break;
			}
			default :{
				XHMM::xhmmOutputManager::printWarning("No centering performed since 'centerType' not given");
				break;
			}
		}
	}

	xhmmOutputManager::printMatrix(*rdMat, outFile);
	xhmmOutputManager::printBasicContainer(*excludedTargets, outputExcludedTargets);
	xhmmOutputManager::printBasicContainer(*excludedSamples, outputExcludedSamples);

	delete rdMat;
	delete excludedTargets;
	delete excludedSamples;
}

HMM_PP::DoubleMat*
XHMM::ProcessMatrixOutputManager::filterTargetProperties(HMM_PP::DoubleMat* rdMat,
		const double minMeanTargetRD, const double maxMeanTargetRD,
		const double minSdTargetRD, const double maxSdTargetRD,
		xhmmInputManager::IntervalSet* excludedTargets) {

	HMM_PP::ullintPair dims = rdMat->getDims();
	const ullint nrow = dims.first;
	const ullint ncol = dims.second;

	set<ullint>* removeColumns = new set<ullint>();

	for (ullint j = 0; j < ncol; ++j) {
		HMM_PP::DistributionStatistics<double> ds;
		for (uint i = 0; i < nrow; ++i)
			ds.observeVal((*rdMat)(i,j));

		double mean = ds.getMean();
		double stdDev = ds.getStdDev();

		if (mean < minMeanTargetRD || mean > maxMeanTargetRD || stdDev < minSdTargetRD || stdDev > maxSdTargetRD) {
			string targStr = rdMat->colName(j);
			cerr << "Removing target " << targStr << " with mean= " << mean << ", sd= " << stdDev << endl;

			removeColumns->insert(j);
			excludedTargets->insert(Interval(targStr));
		}
	}

	if (!removeColumns->empty()) {
		HMM_PP::DoubleMat* newRdMat = rdMat->deleteRowsAndColumns(NULL, removeColumns);
		delete rdMat;
		rdMat = newRdMat;
	}
	delete removeColumns;

	return rdMat;
}

HMM_PP::DoubleMat*
XHMM::ProcessMatrixOutputManager::filterSampleProperties(HMM_PP::DoubleMat* rdMat,
		const double minMeanSampleRD, const double maxMeanSampleRD,
		const double minSdSampleRD, const double maxSdSampleRD,
		xhmmInputManager::StringSet* excludedSamples) {

	HMM_PP::ullintPair dims = rdMat->getDims();
	const ullint nrow = dims.first;
	const ullint ncol = dims.second;

	set<ullint>* removeRows = new set<ullint>();

	for (uint i = 0; i < nrow; ++i) {
		HMM_PP::DistributionStatistics<double> ds;
		for (uint j = 0; j < ncol; ++j)
			ds.observeVal((*rdMat)(i,j));

		double mean = ds.getMean();
		double stdDev = ds.getStdDev();

		if (mean < minMeanSampleRD || mean > maxMeanSampleRD || stdDev < minSdSampleRD || stdDev > maxSdSampleRD) {
			string samp = rdMat->rowName(i);
			cerr << "Removing sample " << samp << " with mean= " << mean << ", sd= " << stdDev << endl;

			removeRows->insert(i);
			excludedSamples->insert(samp);
		}
	}

	if (!removeRows->empty()) {
		HMM_PP::DoubleMat* newRdMat = rdMat->deleteRowsAndColumns(removeRows, NULL);
		delete rdMat;
		rdMat = newRdMat;
	}
	delete removeRows;

	return rdMat;
}
