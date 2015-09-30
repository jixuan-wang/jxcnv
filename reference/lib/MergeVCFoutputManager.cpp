#include "MergeVCFoutputManager.hpp"
#include "xhmmOutputManager.hpp"
#include "AuxiliaryCNVoutputManager.hpp"
#include "xhmmInputManager.hpp"
#include "OnDiskMatrixTransposer.hpp"

#include <Exception.hpp>
#include <NamedMatrix.hpp>

#include <sstream>
#include <list>
#include <string>
using namespace std;

XHMM::MergeVCFoutputManager::MergeVCFoutputManager(list<string>* VCFfiles, list<string>* VCFlistFiles)
: _readers(new list<VCFreader*>()), _samples(new vector<string>()) {

	for (list<string>::const_iterator listIt = VCFlistFiles->begin(); listIt != VCFlistFiles->end(); ++listIt) {
		const string& listFile = *listIt;
		cerr << "Reading list of VCF files from " << listFile << endl;

		list<string>* filesToAdd = xhmmInputManager::readStringsFromFile(listFile);
		VCFfiles->insert(VCFfiles->end(), filesToAdd->begin(), filesToAdd->end());
		delete filesToAdd;
	}
	delete VCFlistFiles;

	if (VCFfiles->empty())
		throw new Exception("Must provide at least one VCF file to merge");

	const VCFheader* firstHeader = NULL;
	set<string>* allSamps = new set<string>();

	for (list<string>::const_iterator fileIt = VCFfiles->begin(); fileIt != VCFfiles->end(); ++fileIt) {
		VCFreader* reader = new VCFreader(*fileIt);

		const VCFheader* curHeader = reader->getHeader();
		if (firstHeader == NULL)
			firstHeader = curHeader;
		else if (*curHeader != *firstHeader)
			throw new Exception("Headers across VCF do not match EXACTLY!");

		_readers->push_back(reader);

		const vector<string>* samps = reader->getSamples();
		for (vector<string>::const_iterator sIt = samps->begin(); sIt != samps->end(); ++sIt) {
			const string& samp = *sIt;
			if (allSamps->find(samp) != allSamps->end())
				throw new Exception("Repeated sample across VCF files being merged: " + samp);
			allSamps->insert(samp);
			_samples->push_back(samp);
		}
	}

	delete VCFfiles;
}

XHMM::MergeVCFoutputManager::~MergeVCFoutputManager() {
	for (list<VCFreader*>::const_iterator it = _readers->begin(); it != _readers->end(); ++it)
		delete *it;
	delete _readers;
}

void XHMM::MergeVCFoutputManager::mergeVCFs(string outFile) {
	if (_readers->empty())
		return;
	const VCFheader* header = (*(_readers->begin()))->getHeader();

	HMM_PP::ostreamWriter* outStream = HMM_PP::utils::getOstreamWriterFromFile(outFile);
	if (outStream == NULL)
		return;
	VCFwriter* out = new VCFwriter(outStream, header, _samples);

	while (true) {
		bool oneHasNext = false;
		bool oneDone = false;
		for (list<VCFreader*>::const_iterator readIt = _readers->begin(); readIt != _readers->end(); ++readIt) {
			VCFreader* reader = *readIt;
			if (reader->hasNext())
				oneHasNext = true;
			else
				oneDone = true;
		}
		if (oneHasNext && oneDone)
			throw new Exception("All VCF files must have the same number of variants");
		if (!oneHasNext)
			break;

		// By necessity, all readers satisfy hasNext():
		Variant* firstVar = NULL;
		VCFwriter::AlleleAndNonAlleleInfo* firstAllInfo = NULL;
		Genotypes* allGenotypes = new Genotypes();

		for (list<VCFreader*>::const_iterator readIt = _readers->begin(); readIt != _readers->end(); ++readIt) {
			VCFreader* reader = *readIt;
			Variant* var = reader->next();
			allGenotypes->insert(var->getGenotypes()->begin(), var->getGenotypes()->end());

			if (firstVar == NULL) {
				firstVar = var;
				firstAllInfo = VCFwriter::splitInfoOnAlleleInfo(firstVar->getInfo()->infoMap());
			}
			else {
				if (var->getChr() != firstVar->getChr()
						|| var->getPos() != firstVar->getPos()
						|| var->getID() != firstVar->getID()
						|| var->getRef() != firstVar->getRef()
						|| var->getAlt() != firstVar->getAlt()
						|| var->getQual() != firstVar->getQual()
						|| var->getFilter() != firstVar->getFilter()
						|| var->getFormat() != firstVar->getFormat())
					throw new Exception("VCF row mismatch");

				VCFwriter::AlleleAndNonAlleleInfo* allInfo = VCFwriter::splitInfoOnAlleleInfo(var->getInfo()->infoMap());
				if (allInfo->info() != firstAllInfo->info())
					throw new Exception("VCF row mismatch in non-allele statistics of INFO field");
				firstAllInfo->allStats() += allInfo->allStats();
				delete allInfo;

				delete var;
			}
		}

		InfoMap* infoMap = new InfoMap(firstAllInfo->info());
		uint an = firstAllInfo->allStats().getAN();
		const list<uint>* ac = firstAllInfo->allStats().getAC();

		list<BaseReal>* af = new list<BaseReal>();
		for (list<uint>::const_iterator acIt = ac->begin(); acIt != ac->end(); ++acIt) {
			BaseReal afVal = 0;
			if (an > 0)
				afVal = *acIt / (BaseReal) an;
			af->push_back(afVal);
		}

		stringstream acStr;
		for (list<uint>::const_iterator acIt = ac->begin(); acIt != ac->end(); ++acIt) {
			if (acIt != ac->begin())
				acStr << ",";
			acStr << *acIt;
		}
		(*infoMap)["AC"] = acStr.str();

		stringstream afStr;
		GenotypeOutputManager::initOstream(afStr);
		for (list<BaseReal>::const_iterator afIt = af->begin(); afIt != af->end(); ++afIt) {
			if (afIt != af->begin())
				afStr << ",";
			afStr << *afIt;
		}
		(*infoMap)["AF"] = afStr.str();

		stringstream anStr;
		anStr << an;
		(*infoMap)["AN"] = anStr.str();

		delete firstAllInfo;

		InfoKeys* infoKeys = new InfoKeys(firstVar->getInfo()->infoKeys());

		Info* info = new Info(infoMap, infoKeys);
		firstVar->setInfo(info);

		firstVar->setGenotypes(allGenotypes);

		*out << *firstVar;
		delete firstVar;
	}

	delete out;
}

XHMM::MergeVCFoutputManager::VCFwriter::VCFwriter(HMM_PP::ostreamWriter* out, const VCFheader* header, vector<string>* samples)
: _out(out), _samples(new vector<string>(*samples)) {

	ostream& stream = (*_out)();
	GenotypeOutputManager::initOstream(stream);

	for (VCFheader::const_iterator headIt = header->begin(); headIt != header->end(); ++headIt)
		stream << *headIt << endl;

	GenotypeOutputManager::printSampleHeaderLine(stream, _samples);
}

XHMM::MergeVCFoutputManager::VCFwriter::~VCFwriter() {
	delete _out;
	delete _samples;
}

XHMM::MergeVCFoutputManager::VCFwriter& XHMM::MergeVCFoutputManager::VCFwriter::operator<<(Variant& v) {
	ostream& stream = (*_out)();

	stream
	<< v.getChr()
	<< '\t' << v.getPos()
	<< '\t' << v.getID()
	<< '\t' << v.getRef()
	<< '\t' << v.getAlt()
	<< '\t' << v.getQual()
	<< '\t' << v.getFilter();

	stream
	<< '\t';
	const Info* info = v.getInfo();
	const InfoMap& infoMap = info->infoMap();
	const InfoKeys& infoKeys = info->infoKeys();
	for (InfoKeys::const_iterator inIt = infoKeys.begin(); inIt != infoKeys.end(); ++inIt) {
		const string& key = *inIt;
		const string& val = infoMap.find(key)->second;

		if (inIt != infoKeys.begin())
			stream << ';';
		stream << key;
		if (val != "")
			stream << '=' << val;
	}

	stream
	<< '\t' << v.getFormat();

	const Genotypes* gt = v.getGenotypes();
	for (vector<string>::const_iterator sampIt = _samples->begin(); sampIt != _samples->end(); ++sampIt) {
		string gtStr = XHMM::GenotypeOutputManager::MISSING_VAL;
		Genotypes::const_iterator findSamp = gt->find(*sampIt);
		if (findSamp != gt->end())
			gtStr = findSamp->second;
		stream << '\t' << gtStr;
	}

	stream << endl;

	return *this;
}

XHMM::MergeVCFoutputManager::VCFreader::VCFreader(string VCFfile)
: _stream(HMM_PP::utils::getIstreamLineReaderFromFile(VCFfile)),
  _header(new VCFheader()), _samples(new vector<string>()) {
	readHeaderAndSamples();
}

XHMM::MergeVCFoutputManager::VCFreader::~VCFreader() {
	delete _stream;
	delete _header;
	delete _samples;
}

void XHMM::MergeVCFoutputManager::VCFreader::readHeaderAndSamples() {
	while (!_stream->eof()) {
		string* nextRow = new string();
		*_stream >> *nextRow;

		if (nextRow->substr(0, GenotypeOutputManager::HEADER_ROW_PREFIX.size()) == GenotypeOutputManager::HEADER_ROW_PREFIX) {
			_header->push_back(*nextRow);
			delete nextRow;
		}
		else {
			stringstream* columnsRowStream = new stringstream(*nextRow);

			vector<string>* fields = new vector<string>();
			while (*columnsRowStream && !columnsRowStream->eof()) {
				string field;
				// Instead of splitting by whitespace, use ONLY tab as delimiter for the header line (to allow for sample names with space in them):
				if (!getline(*columnsRowStream, field, '\t'))
					break;

				fields->push_back(field);
			}
			delete columnsRowStream;

			uint ind = 0;
			for (list<string>::const_iterator i = GenotypeOutputManager::VCF_COLUMNS.begin(); i != GenotypeOutputManager::VCF_COLUMNS.end(); ++i) {
				if (*i != (*fields)[ind++])
					throw new Exception("Malformed VCF header: " + *nextRow);
			}
			delete nextRow;

			// Get samples:
			for (; ind < fields->size(); ++ind)
				_samples->push_back((*fields)[ind]);
			delete fields;

			break;
		}
	}
}

XHMM::MergeVCFoutputManager::Variant::Variant(const string& variantRow, const vector<string>* samples)
: _vals(new vector<string>()), _samples(samples),
  _pos(-1), _info(NULL), _genotypes(NULL) {
	stringstream* columnVals = new stringstream(variantRow);
	while (*columnVals && !columnVals->eof()) {
		string val;
		// Instead of splitting by whitespace, use ONLY tab as delimiter for the header line (to allow for sample names with space in them):
		if (!getline(*columnVals, val, '\t'))
			break;

		_vals->push_back(val);
	}
	delete columnVals;

	if (_vals->size() != GenotypeOutputManager::VCF_COLUMNS.size() + _samples->size())
		throw new Exception("Invalid # of columns in VCF row: " + variantRow);
}

XHMM::MergeVCFoutputManager::Variant::~Variant() {
	delete _vals;

	if (_info != NULL)
		delete _info;

	if (_genotypes != NULL)
		delete _genotypes;
}

int XHMM::MergeVCFoutputManager::Variant::getPos() {
	if (_pos < 0) {
		stringstream str((*_vals)[GenotypeOutputManager::POS]);
		str >> _pos;
		if (!str)
			throw new Exception("Invalid VCF position " + (*_vals)[GenotypeOutputManager::POS]);
	}
	return _pos;
}

const XHMM::MergeVCFoutputManager::Info* XHMM::MergeVCFoutputManager::Variant::getInfo() {
	if (_info == NULL) {
		InfoMap* infoMap = new InfoMap();
		InfoKeys* infoKeys = new InfoKeys();

		stringstream* infoKeyVals = new stringstream((*_vals)[GenotypeOutputManager::INFO]);
		while (*infoKeyVals && !infoKeyVals->eof()) {
			string keyValPair;
			if (!getline(*infoKeyVals, keyValPair, ';'))
				break;

			string key;
			string val;

			stringstream keyValPairStr(keyValPair);
			if (!getline(keyValPairStr, key, '='))
				throw new Exception("Invalid key=value pair in INFO field: " + keyValPair);
			keyValPairStr >> val;
			if (!keyValPairStr)
				val = "";

			(*infoMap)[key] = val;
			infoKeys->push_back(key);
		}
		delete infoKeyVals;

		_info = new Info(infoMap, infoKeys);
	}
	return _info;
}

XHMM::MergeVCFoutputManager::VCFwriter::AlleleStats::AlleleStats(const string& ac, const string& an)
: _ac(new list<uint>()), _an(0) {
	stringstream* acVals = new stringstream(ac);
	while (*acVals && !acVals->eof()) {
		string acValString;
		if (!getline(*acVals, acValString, ','))
			break;
		stringstream acValStream(acValString);

		uint acVal;
		acValStream >> acVal;
		if (!acValStream)
			throw new Exception("Unable to read numeric value from AC vector");
		_ac->push_back(acVal);
	}
	delete acVals;

	stringstream anStream(an);
	anStream >> _an;
	if (!anStream)
		throw new Exception("Unable to read numeric value for AN");
}

XHMM::MergeVCFoutputManager::VCFwriter::AlleleStats::~AlleleStats() {
	delete _ac;
}

XHMM::MergeVCFoutputManager::VCFwriter::AlleleStats& XHMM::MergeVCFoutputManager::VCFwriter::AlleleStats::operator+=(const AlleleStats& add) {
	if (add._ac->size() != _ac->size())
		throw new Exception("Cannot add AC vectors of different lengths");

	list<uint>::iterator thisIt = _ac->begin();
	list<uint>::iterator addIt = add._ac->begin();
	while (true) {
		if (thisIt == _ac->end())
			break;

		*thisIt += *addIt;

		++thisIt;
		++addIt;
	}

	_an += add._an;

	return *this;
}

XHMM::MergeVCFoutputManager::VCFwriter::AlleleAndNonAlleleInfo* XHMM::MergeVCFoutputManager::VCFwriter::splitInfoOnAlleleInfo(const InfoMap& info) {
	string ac;
	string an;

	InfoMap* otherInfo = new InfoMap(info);

	InfoMap::iterator acIt = otherInfo->find("AC");
	if (acIt != otherInfo->end()) {
		ac = acIt->second;
		otherInfo->erase(acIt);
	}
	else
		throw new Exception("Missing AC in INFO");

	InfoMap::iterator afIt = otherInfo->find("AF");
	if (afIt != otherInfo->end())
		otherInfo->erase(afIt);
	else
		throw new Exception("Missing AF in INFO");

	InfoMap::iterator anIt = otherInfo->find("AN");
	if (anIt != otherInfo->end()) {
		an = anIt->second;
		otherInfo->erase(anIt);
	}
	else
		throw new Exception("Missing AN in INFO");

	AlleleStats* alleleStats = new AlleleStats(ac, an);

	return new AlleleAndNonAlleleInfo(alleleStats, otherInfo);
}

const XHMM::MergeVCFoutputManager::Genotypes* XHMM::MergeVCFoutputManager::Variant::getGenotypes() {
	if (_genotypes == NULL) {
		_genotypes = new Genotypes();

		for (uint i = 0; i < _samples->size(); ++i)
			(*_genotypes)[(*_samples)[i]] = (*_vals)[GenotypeOutputManager::VCF_COLUMNS.size() + i];
	}
	return _genotypes;
}

XHMM::MergeVCFoutputManager::Variant* XHMM::MergeVCFoutputManager::VCFreader::next() {
	if (!hasNext())
		throw new Exception("Cannot read additional variants");

	string* nextRow = new string();
	*_stream >> *nextRow;
	Variant* v = new Variant(*nextRow, _samples);
	delete nextRow;

	return v;
}
