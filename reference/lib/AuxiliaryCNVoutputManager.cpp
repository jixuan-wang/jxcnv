#include "AuxiliaryCNVoutputManager.hpp"
#include "SampleGenotypeQualsCalculator.hpp"
#include "CNVmodelParams.hpp"
#include "Genotype.hpp"
#include "AlleleQuals.hpp"
#include "DiscoverOutputManager.hpp"

#include <params.hpp>
#include <Model.hpp>
#include <Data.hpp>
#include <InferredDataFit.hpp>

#include <set>
#include <fstream>
using namespace std;

XHMM::AuxiliaryCNVoutputManager::AuxiliaryCNVoutputManager
(string auxOutputFile,
		uint auxUpstreamPrintTargs, uint auxDownstreamPrintTargs,
		bool printOrigData)
: _auxStream(HMM_PP::utils::getOstreamWriterFromFile(auxOutputFile)),
  _auxUpstreamPrintTargs(auxUpstreamPrintTargs), _auxDownstreamPrintTargs(auxDownstreamPrintTargs),
  _model(NULL), _dataLoader(NULL), _printOrigData(printOrigData) {

	if (_auxStream == NULL)
		throw new Exception("Cannot give empty file name for auxiliary output");

	(*_auxStream)() << setiosflags(ios::fixed);
	(*_auxStream)() << setprecision(POSTERIOR_AND_RD_PRECISION);

	printTargetHeader();
}

XHMM::AuxiliaryCNVoutputManager::~AuxiliaryCNVoutputManager() {
	delete _auxStream;
}

void XHMM::AuxiliaryCNVoutputManager::setModelAndReadDepthLoader(HMM_PP::Model* model, const ReadDepthMatrixLoader* dataLoader) {
	if (_model != NULL || _dataLoader != NULL)
		throw new Exception("Cannot call setModelAndReadDepthLoader() more than once");
	_model = model;
	_dataLoader = dataLoader;
}

void XHMM::AuxiliaryCNVoutputManager::printTargetHeader() {
	ostream& outStream = (*_auxStream)();

	outStream
	<< "SAMPLE"
	<< '\t' << XHMM::DiscoverOutputManager::CNV
	<< '\t' << "FULL_" << XHMM::DiscoverOutputManager::INTERVAL
	<< '\t' << "TARGET_IND"
	<< '\t' << "TARGET"
	<< '\t' << "CHR"
	<< '\t' << "MID_BP"
	<< '\t' << "POSTERIOR"
	<< '\t' << "RD";

	if (_printOrigData)
		outStream
		<< '\t' << "ORIG_RD";

	outStream
	<< endl;
}

void XHMM::AuxiliaryCNVoutputManager::printCNVtargets(SampleGenotypeQualsCalculator* genotyper, const uint t1, const uint t2, const vector<uint>& types, HMM_PP::Data* origData) {
	if (_model == NULL || _dataLoader == NULL)
		throw new Exception("Must call setModelAndReadDepthLoader() before using object");

	const HMM_PP::InferredDataFit* dataFit = genotyper->getDataFit();

	const Interval cnvTarg = genotyper->getMergedInterval(t1, t2);

	const uint chrStart = _dataLoader->getChrStartTargetIndex(cnvTarg.getChr());
	const uint chrStop = _dataLoader->getChrStopTargetIndex(cnvTarg.getChr());

	uint start;
	if (t1 >= chrStart + _auxUpstreamPrintTargs)
		start = t1 - _auxUpstreamPrintTargs;
	else
		start = chrStart;

	uint stop;
	if (t2 + _auxDownstreamPrintTargs <= chrStop)
		stop = t2 + _auxDownstreamPrintTargs;
	else
		stop = chrStop;

	ostream& outStream = (*_auxStream)();

	for (uint t = start; t <= stop; ++t) {
		const Interval& targT = _dataLoader->getTarget(t);
		HMM_PP::ModelParams::DataVal* rd = genotyper->calcMeanRD(t, t);

		stringstream targInd;
		if (t < t1)
			targInd << "U-" << (t1 - t);
		else if (t > t2)
			targInd << "D+" << (t - t2);
		else
			targInd << (t - t1 + 1);

		outStream
		<< dataFit->getData()->getId()
		<< '\t';

		const HMM_PP::ModelParams* params = dataFit->getModel()->getModelParams();
		for (uint i = 0; i < types.size(); ++i) {
			if (i > 0)
				outStream << ",";
			outStream
			<< params->state(types[i]);
		}

		outStream
		<< '\t' << cnvTarg
		<< '\t' << targInd.str()
		<< '\t' << targT
		<< '\t' << targT.getChr()
		<< '\t' << targT.midpoint()
		<< '\t';

		for (uint i = 0; i < types.size(); ++i) {
			real posterior = dataFit->getGamma(t, types[i]); // == data->calcPosterior(t, t, types[i])
			if (i > 0)
				outStream << ",";
			outStream
			<< posterior;
		}

		outStream
		<< '\t' << *rd;

		delete rd;

		if (_printOrigData) {
			outStream
			<< '\t';

			if (origData != NULL) {
				HMM_PP::ModelParams::DataVal* origRD = _model->calcRepresentativeDataVal(origData, t, t);
				outStream
				<< *origRD;
				delete origRD;
			}
			else
				outStream
				<< "NA";
		}

		outStream
		<< endl;
	}
}

void XHMM::AuxiliaryCNVoutputManager::printCNVtargets(SampleGenotypeQualsCalculator* genotyper, const uint t1, const uint t2, const uint type, HMM_PP::Data* origData) {
	vector<uint> types;
	types.push_back(type);

	printCNVtargets(genotyper, t1, t2, types, origData);
}
