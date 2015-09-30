#include "CNVmodelParams.hpp"
#include "xhmmDriver.hpp"
#include "xhmmInputManager.hpp"
#include "xhmmOutputManager.hpp"
#include "DiscoverOutputManager.hpp"
#include "ReadDepthMatrixLoader.hpp"
#include "Interval.hpp"
#include "AuxiliaryCNVoutputManager.hpp"
#include "ProcessMatrixOutputManager.hpp"
#include "PCA_NormalizeOutputManager.hpp"
#include "PrepareTargetsManager.hpp"
#include "GATKdepthMergerOutputManager.hpp"
#include "GenotypeOutputManager.hpp"
#include "MatrixDecomp.hpp"
#include "MergeVCFoutputManager.hpp"

#include <utils.hpp>
#include <Model.hpp>
#include <Exception.hpp>
#include <ModelWithDB.hpp>
#include <ModelFromStream.hpp>
#include <ModelFromStreamNoStorage.hpp>
#include <PreciseNonNegativeReal.hpp>

#include <libgen.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
using namespace std;

static const double D_INF = HMM_PP::PreciseNonNegativeReal<double>::REAL_INFINITY;

XHMM::xhmmDriver::xhmmDriver(const gengetopt_args_info& args_info)
: _args_info(args_info),
  _paramStream(NULL), _dbFile(NULL),
  _model(NULL), _dataLoader(NULL),
  _origDataLoader(NULL), _auxOutput(NULL),
  _keepIDs(NULL) {}

XHMM::xhmmDriver::~xhmmDriver() {
	if (_paramStream != NULL)
		delete _paramStream;

	if (_dbFile != NULL)
		delete _dbFile;

	if (_model != NULL)
		delete _model;

	if (_dataLoader != NULL)
		delete _dataLoader;

	if (_origDataLoader != NULL)
		delete _origDataLoader;

	if (_auxOutput != NULL)
		delete _auxOutput;

	if (_keepIDs != NULL)
		delete _keepIDs;
}

bool XHMM::xhmmDriver::run() {
	runCommonInits();

	if (_args_info.prepareTargets_given)
		return prepareTargets();

	if (_args_info.mergeGATKdepths_given)
		return mergeGATKdepths();

	if (_args_info.matrix_given)
		return processMatrix();

	if (_args_info.PCA_given)
		return PCA();

	if (_args_info.normalize_given)
		return normalize();

	if (_args_info.discover_given)
		return discover();

	if (_args_info.genotype_given)
		return genotype();

	if (_args_info.mergeVCFs_given)
		return mergeVCFs();

	if (_args_info.transition_given)
		return transition();

	if (_args_info.createDB_given)
		return createDB();

	if (_args_info.printHMM_given)
		return printHMM();

	throw new Exception("No mode parameter given");

	return false;
}

void XHMM::xhmmDriver::runCommonInits() {
	HMM_PP::MatrixDecomp::initSignalHandlers();

	setChrOrderFromReference();

	readSampleIDsToKeep();
}

void XHMM::xhmmDriver::setChrOrderFromReference() {
	if (_args_info.referenceFASTA_given) {
		xhmmInputManager::Table* fastaIndexTable = xhmmInputManager::readFastaIndexTable(_args_info.referenceFASTA_arg);

		XHMM::Interval::CHR_TO_INDEX.clear();
		for (uint i = 0; i < fastaIndexTable->getNumRows(); ++i) // Chromosomes are in the 0-th column:
			XHMM::Interval::CHR_TO_INDEX[fastaIndexTable->getEntry(i, 0)] = i;

		delete fastaIndexTable;
	}
}

void XHMM::xhmmDriver::readSampleIDsToKeep() {
	if (_args_info.keepSampleIDs_given) {
		xhmmInputManager::Table* idsTable = xhmmInputManager::readTable(_args_info.keepSampleIDs_arg, false);
		_keepIDs = xhmmInputManager::tableColumnToSet(idsTable, 1);
	}
}

bool XHMM::xhmmDriver::prepareTargets() {
	PrepareTargetsManager* ptm = new PrepareTargetsManager();

	ptm->prepareTargets(xhmmInputManager::getListFromArray(_args_info.targets_arg, _args_info.targets_given), _args_info.mergedTargets_arg);

	delete ptm;

	return true;
}

bool XHMM::xhmmDriver::mergeGATKdepths() {
	xhmmInputManager::OneToOneStrMap* sampMap = NULL;
	if (_args_info.sampleIDmap_given)
		sampMap = xhmmInputManager::tableToMap(xhmmInputManager::readTable(_args_info.sampleIDmap_arg, false), _args_info.fromID_arg, _args_info.toID_arg);

	GATKdepthMergerOutputManager* gdm =
			new GATKdepthMergerOutputManager(xhmmInputManager::getListFromArray(_args_info.GATKdepths_arg, _args_info.GATKdepths_given),
					xhmmInputManager::getListFromArray(_args_info.GATKdepthsList_arg, _args_info.GATKdepthsList_given), sampMap, _args_info.columnSuffix_arg);

	gdm->mergeGATKdepths(_args_info.outputMatrix_arg, !_args_info.outputTargetsBySamples_flag, _args_info.rdPrecision_arg);

	delete gdm;

	return true;
}

bool XHMM::xhmmDriver::processMatrix() {
	ProcessMatrixOutputManager* pm = new ProcessMatrixOutputManager();

	if (_args_info.log10_given && _args_info.log10_arg <= 0)
		throw new Exception("Must provide log10 argument value > 0");

	pm->processReadDepthMatrix(_args_info.readDepths_arg, _args_info.outputMatrix_arg,
			xhmmInputManager::getListFromArray(_args_info.excludeTargets_arg, _args_info.excludeTargets_given),
			xhmmInputManager::getListFromArray(_args_info.excludeChromosomeTargets_arg, _args_info.excludeChromosomeTargets_given),
			xhmmInputManager::getListFromArray(_args_info.excludeSamples_arg, _args_info.excludeSamples_given),
			_args_info.minTargetSize_arg, (_args_info.maxTargetSize_given ? _args_info.maxTargetSize_arg : UINT_INFINITY),
			(_args_info.minMeanTargetRD_given ? _args_info.minMeanTargetRD_arg : - D_INF), (_args_info.maxMeanTargetRD_given ? _args_info.maxMeanTargetRD_arg : D_INF),
			_args_info.minSdTargetRD_arg, (_args_info.maxSdTargetRD_given ? _args_info.maxSdTargetRD_arg : D_INF),
			(_args_info.minMeanSampleRD_given ? _args_info.minMeanSampleRD_arg: - D_INF), (_args_info.maxMeanSampleRD_given ? _args_info.maxMeanSampleRD_arg : D_INF),
			_args_info.minSdSampleRD_arg, (_args_info.maxSdSampleRD_given ? _args_info.maxSdSampleRD_arg : D_INF),
			_args_info.scaleDataBySum_flag, _args_info.scaleDataBySumType_arg, _args_info.scaleDataBySumFactor_arg,
			_args_info.log10_given, (_args_info.log10_given ? _args_info.log10_arg : 0),
			_args_info.centerData_flag, _args_info.centerType_arg,
			_args_info.zScoreData_flag,
			(_args_info.outputExcludedTargets_given ? _args_info.outputExcludedTargets_arg : ""),
			(_args_info.outputExcludedSamples_given ? _args_info.outputExcludedSamples_arg : ""));

	delete pm;

	return true;
}

bool XHMM::xhmmDriver::PCA() {
	string workDir = "";
	if (_args_info.PCA_saveMemory_given) {
		workDir = _args_info.PCA_saveMemory_arg;
		if (workDir == "") {
			char* outFile = strdup(_args_info.PCAfiles_arg);
			workDir = dirname(outFile);
			free(outFile);
		}
	}

	PCA_NormalizeOutputManager* pcaManager = new PCA_NormalizeOutputManager(_args_info.PCAfiles_arg, workDir);
	pcaManager->PCA(_args_info.readDepths_arg);
	delete pcaManager;

	return true;
}

bool XHMM::xhmmDriver::normalize() {
	PCA_NormalizeOutputManager* pcaManager = new PCA_NormalizeOutputManager(_args_info.PCAfiles_arg);
	pcaManager->normalize(_args_info.readDepths_arg, _args_info.normalizeOutput_arg, _args_info.PCnormalizeMethod_arg, _args_info.numPCtoRemove_arg, _args_info.PVE_mean_factor_arg, _args_info.PVE_contrib_arg);
	delete pcaManager;

	return true;
}

bool XHMM::xhmmDriver::discover() {
	createOrigDataLoaderAndAuxOutput();

	string xcnvFile = _args_info.xcnv_arg;

	string posteriorFiles = "";
	if (_args_info.posteriorFiles_given)
		posteriorFiles = _args_info.posteriorFiles_arg;
	else if (xcnvFile != HMM_PP::utils::STD_STREAM)
		posteriorFiles = xcnvFile;

	DiscoverOutputManager* outputManager = new DiscoverOutputManager(xcnvFile, _auxOutput, posteriorFiles, _args_info.discoverSomeQualThresh_arg, _args_info.maxQualScore_arg, _args_info.scorePrecision_arg, _origDataLoader);
	createModelAndLoader(outputManager);
	outputManager->setModelAndReadDepthLoader(_model, _dataLoader);

	// First, try to estimate HMM parameters:
	if (_args_info.optDiscover_given)
		_model->nelderMead();

	// Fit data to model with fixed parameters:
	const bool CALC_DIGAMMAS = false; // since we don't need all of them, and the quality calculator will calculate what is needed on the fly for the start and stop quals
	const bool RUN_VITERBI = true; // since we want the max-probability sequences
	real lik = _model->fit(CALC_DIGAMMAS, RUN_VITERBI);

	cerr
	<< endl
	<< "log10-likelihood of all data: " << lik.getLog10Value() << endl;

	// Output results:
#if defined(DEBUG)
	_model->displayFittedSequence(cerr);
#endif
	// NOTE: this command does nothing if output was already done on-the-fly by ModelFromStreamNoStorage
	outputManager->printAllCNVs(_model);

	delete outputManager;

	return true;
}

bool XHMM::xhmmDriver::genotype() {
	const bool genotypeSubsegments = _args_info.subsegments_flag;
	const uint maxTargetsInSubsegment = _args_info.maxTargetsInSubsegment_arg;

	xhmmInputManager::Table* cnvTable = xhmmInputManager::readTable(_args_info.gxcnv_arg);
	cerr << "Read " << cnvTable->getNumRows() << " CNV lines from " << _args_info.gxcnv_arg << endl;

	const DiscoverOutputManager::IntervalToThreshMap* intervals = DiscoverOutputManager::CNVtableToIntervals(cnvTable);
	cerr << "Will genotype " << intervals->size() << " unique intervals";
	if (genotypeSubsegments) cerr << " (and their subsegments overlapping up to " << maxTargetsInSubsegment << " targets)";
	cerr << endl;

	createOrigDataLoaderAndAuxOutput();

	GenotypeOutputManager* outputManager = new GenotypeOutputManager(_args_info.vcf_arg, intervals, _auxOutput, _args_info.genotypeQualThresholdWhenNoExact_arg, genotypeSubsegments, maxTargetsInSubsegment, _origDataLoader, _args_info.maxQualScore_arg, _args_info.scorePrecision_arg);
	createModelAndLoader(outputManager);
	outputManager->setModelAndReadDepthLoader(_model, _dataLoader);

	// Fit data to model with fixed parameters, and use the results to genotype the CNVs in all the samples:
	_model->fit(false, false, true); // since we don't need all digammas [the quality calculator will calculate what is needed on the fly for the start and stop quals] and we don't need Viterbi, but we do want calculations to be timed

	delete outputManager;

	return true;
}

bool XHMM::xhmmDriver::mergeVCFs() {
	MergeVCFoutputManager* mvm =
			new MergeVCFoutputManager(xhmmInputManager::getListFromArray(_args_info.mergeVCF_arg, _args_info.mergeVCF_given),
					xhmmInputManager::getListFromArray(_args_info.mergeVCFlist_arg, _args_info.mergeVCFlist_given));

	mvm->mergeVCFs(_args_info.vcf_arg);

	delete mvm;

	return true;
}

bool XHMM::xhmmDriver::printHMM() {
	CNVmodelParams* cnvParams = new CNVmodelParams(*openParamFile());
	xhmmOutputManager::printInitProbs(cnvParams, cerr);
	delete cnvParams;

	return false;
}

bool XHMM::xhmmDriver::transition() {
	CNVmodelParams* cnvParams = new CNVmodelParams(*openParamFile());

	xhmmOutputManager::printInitProbs(cnvParams, cerr);

	const string distFile = HMM_PP::utils::STD_STREAM;
	xhmmOutputManager::printTransitionMatricesForDistancesFile(cnvParams, distFile, cerr);

	delete cnvParams;

	return false;
}

bool XHMM::xhmmDriver::createDB() {
	if (getDBfile() == NULL)
		throw new Exception("Must provide database file in which to load data");

	createModelAndLoader();

	_model->loadAllData();

	return true;
}

istream* XHMM::xhmmDriver::openParamFile() {
	if (_paramStream == NULL) {
		if (!_args_info.paramFile_given)
			throw new Exception("paramFile option not given");
		string paramFile = _args_info.paramFile_arg;
		_paramStream = new ifstream(paramFile.c_str());
		if (!*_paramStream)
			throw new Exception("Unable to open parameter file: " + paramFile);
	}

	return _paramStream;
}

string* XHMM::xhmmDriver::getDBfile() {
	if (_dbFile == NULL) {
		if (_args_info.DB_given)
			_dbFile = new string(_args_info.DB_arg);
	}

	return _dbFile;
}

void XHMM::xhmmDriver::createModelAndLoader(HMM_PP::DataOutputter<HMM_PP::InferredDataFit>* dataOuputter) {
	// Use the non-homogeneous CNV HMM:
	CNVmodelParams* cnvParams = NULL;
	if (getDBfile() == NULL) {
		string readDepthFile = _args_info.readDepths_arg;
		double maxNormalizedReadDepthVal = _args_info.maxNormalizedReadDepthVal_arg;
		_dataLoader = new ReadDepthMatrixLoader(readDepthFile, _keepIDs, maxNormalizedReadDepthVal);

		cnvParams = new CNVmodelParams(*openParamFile(), _dataLoader);
	}
	else {
		/* TODO: Currently, the targets are not stored in the DB, but they need to be retrieved from SOMEWHERE.
		 * This is complicated by the fact that ModelWithDB only knows about Data and not the meta-data (e.g., target information).
		 */
		throw new Exception("STILL need to create a _dataLoader object that has only targets in it: Retrieve this from DB???");
	}

	if (getDBfile() == NULL) {
		if (((_args_info.discover_given && !_args_info.optDiscover_given) || _args_info.genotype_given) && dataOuputter != NULL) // no need to store data in memory:
			_model = new HMM_PP::ModelFromStreamNoStorage(cnvParams, dataOuputter);
		else
			_model = new HMM_PP::ModelFromStream(cnvParams);
	}
	else
		_model = new HMM_PP::ModelWithDB(cnvParams, *getDBfile());

#if defined(DEBUG)
	_model->printModel(cerr) << endl;
	_model->printExplicitModel(cerr, _dataLoader->getNumTargets() - 1);
#endif

	if (_auxOutput != NULL)
		_auxOutput->setModelAndReadDepthLoader(_model, _dataLoader);
}

void XHMM::xhmmDriver::createOrigDataLoaderAndAuxOutput() {
	// NOTE: original read depths are NOT capped by maxReadDepthVal:
	if (_args_info.origReadDepths_given)
		_origDataLoader = new ReadDepthMatrixLoader(_args_info.origReadDepths_arg, _keepIDs);

	if (_args_info.aux_xcnv_given) {
		const int auxUpstreamPrintTargs = _args_info.auxUpstreamPrintTargs_arg;
		const int auxDownstreamPrintTargs = _args_info.auxDownstreamPrintTargs_arg;

		if (auxUpstreamPrintTargs < 0)
			throw new Exception("auxUpstreamPrintTargs < 0");
		if (auxDownstreamPrintTargs < 0)
			throw new Exception("auxDownstreamPrintTargs < 0");
		_auxOutput = new AuxiliaryCNVoutputManager(_args_info.aux_xcnv_arg, auxUpstreamPrintTargs, auxDownstreamPrintTargs, _origDataLoader != NULL);
	}
}
