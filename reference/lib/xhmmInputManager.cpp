#include "xhmmInputManager.hpp"
#include "utils.hpp"
#include "ReadDepthMatrixLoader.hpp"
#include "xhmmOutputManager.hpp"

#include <Exception.hpp>
#include <NamedMatrix.hpp>
#include <DataLoader.hpp>

#include <sstream>
using namespace std;

const string XHMM::xhmmInputManager::FASTA_SUFFIX = ".fasta";
const string XHMM::xhmmInputManager::FASTA_INDEX_SUFFIX = ".fai";

set<char> XHMM::xhmmInputManager::INTERVALS_EXCLUDE_LINE_CHARS = createIntervalsExcludeCharSet();

set<char> XHMM::xhmmInputManager::createIntervalsExcludeCharSet() {
	set<char> exclude;
	exclude.insert(COMMENT_LINE_CHAR);
	exclude.insert(INTERVALS_HEADER_LINE_CHAR);

	return exclude;
}

XHMM::xhmmInputManager::Table::Table(Row* colNames)
: _colNames(colNames), _colNameToIndex(new map<string, ullint>()), _rows(new vector<Row*>()) {

	for (ullint i = 0; i < colNames->size(); ++i)
		(*_colNameToIndex)[(*_colNames)[i]] = i;
}

XHMM::xhmmInputManager::Table::Table()
: _colNames(new Row()), _colNameToIndex(new map<string, ullint>()), _rows(new vector<Row*>()) {
}

XHMM::xhmmInputManager::Table::~Table() {
	delete _colNames;
	delete _colNameToIndex;
	for (vector<Row*>::const_iterator it = _rows->begin(); it != _rows->end(); ++it)
		delete *it;
	delete _rows;
}

const string& XHMM::xhmmInputManager::Table::getEntry(ullint row, const string& column) const {
	if (row > _rows->size()) {
		stringstream str;
		str << "Invalid table row: " << row;
		throw new Exception(str.str());
	}
	if (!hasColumn(column))
		throw new Exception("Invalid table column: " + column);

	return (*(*_rows)[row])[(*_colNameToIndex)[column]];
}

const string& XHMM::xhmmInputManager::Table::getEntry(ullint row, ullint column) const {
	if (row > _rows->size()) {
		stringstream str;
		str << "Invalid table row: " << row;
		throw new Exception(str.str());
	}
	if (column > _colNames->size()) {
		stringstream str;
		str << "Invalid table column: " << column;
		throw new Exception(str.str());
	}

	return (*(*_rows)[row])[column];
}

void XHMM::xhmmInputManager::Table::addRow(Row* row) {
	if (_colNames->empty()) { // for constructor of Table with NO header
		for (ullint i = 0; i < row->size(); ++i) {
			stringstream str;
			str << i;
			string iString = str.str();
			_colNames->push_back(iString);
			(*_colNameToIndex)[iString] = i;
		}
	}

	if (row->size() != _colNames->size()) {
		stringstream str;
		str << "Input rows must be of equal size: " << _colNames->size();
		throw new Exception(str.str());
	}

	_rows->push_back(row);
}

XHMM::xhmmInputManager::Table* XHMM::xhmmInputManager::readTable(string tableFile, bool header) {
#if defined(DEBUG)
	cerr << "Reading file " << tableFile << endl << endl;
#endif

	HMM_PP::istreamLineReader* tableStream = HMM_PP::utils::getIstreamLineReaderFromFile(tableFile);
	if (tableStream == NULL || tableStream->eof())
		throw new Exception("Unable to read from (possibly empty) file '" + tableFile + "'");

	Table* table = NULL;
	if (header) {
		Table::Row* columns = readRowFromStream(tableStream);
		table = new Table(columns);
	}
	else
		table = new Table();

	while (!tableStream->eof()) {
		Table::Row* row = readRowFromStream(tableStream);
		table->addRow(row);
	}

	delete tableStream;

	return table;
}

XHMM::xhmmInputManager::Table::Row* XHMM::xhmmInputManager::readRowFromStream(HMM_PP::istreamLineReader* tableStream) {
	if (tableStream->eof())
		return NULL;

	string* line = new string();
	*tableStream >> *line;
	stringstream* lineStream = new stringstream(*line);
	delete line;

	Table::Row* row = new Table::Row();
	while (*lineStream && !lineStream->eof()) {
		string val;

		// Instead of splitting by whitespace, use ONLY tab as delimiter for the header line (to allow for sample names with space in them):
		if (!getline(*lineStream, val, '\t'))
			break;

		row->push_back(val);
	}
	delete lineStream;

	return row;
}

XHMM::xhmmInputManager::Table* XHMM::xhmmInputManager::readFastaIndexTable(string fastaFile) {
	string fastaIndexFile = fastaFile + FASTA_INDEX_SUFFIX;

	Table* fastaIndexTable = NULL;
	try {
		fastaIndexTable = readTable(fastaIndexFile, false);
	}
	catch (FileNotFoundException* e) {
		string::size_type expectedSuffixLoc = fastaFile.size() - FASTA_SUFFIX.size();
		if (fastaFile.find(FASTA_SUFFIX) == expectedSuffixLoc) {
			xhmmOutputManager::printWarning(e->getMessage());
			fastaIndexFile = fastaFile.substr(0, expectedSuffixLoc) + FASTA_INDEX_SUFFIX;
			cerr << "Searching for index file " << fastaIndexFile << endl << endl;
			fastaIndexTable = readTable(fastaIndexFile, false);
		}
		else
			throw e;
	}
	catch (Exception* e) {
		stringstream str;
		str << "In reading " << fastaFile << ": " << e->getMessage();
		Exception* throwE = new Exception(str.str());
		delete e;
		throw throwE;
	}

	return fastaIndexTable;
}

list<string>* XHMM::xhmmInputManager::getListFromArray(char** vals, const ullint length) {
	list<string>* valsList = new list<string>();

	for (ullint i = 0; i < length; ++i)
		valsList->push_back(vals[i]);

	return valsList;
}

XHMM::xhmmInputManager::LoadedReadDepths
XHMM::xhmmInputManager::loadReadDepthsFromFile(string readDepthFile, IntervalSet* excludeTargets, StringSet* excludeTargetChromosomes, StringSet* excludeSamples,
		const ullint minTargetSize, const ullint maxTargetSize) {

	HMM_PP::DoubleMat* rdMat = MatrixReader<double>::readMatrixFromFile(readDepthFile);
	StringSet* excludedSamples = new StringSet();
	IntervalSet* excludedTargets = new IntervalSet();

	set<ullint>* excludeSampleIndices = new set<ullint>();
	if (excludeSamples != NULL) {
		for (ullint row = 0; row < rdMat->nrow(); ++row) {
			string samp = rdMat->rowName(row);

			if (excludeSamples->find(samp) != excludeSamples->end()) {
				cerr << "Excluded sample " << samp << endl;

				excludedSamples->insert(samp);
				excludeSampleIndices->insert(row);
			}
		}
	}

	set<ullint>* excludeTargetIndices = new set<ullint>();
	if (excludeTargets != NULL || excludeTargetChromosomes != NULL || minTargetSize > 0 || maxTargetSize < ULLINT_INFINITY) {
		for (ullint j = 0; j < rdMat->ncol(); ++j) {
			const Interval targJ(rdMat->colName(j));
			const ullint targLen = targJ.span();
			bool targLenFails = (targLen < minTargetSize || targLen > maxTargetSize);

			if ((excludeTargets != NULL && excludeTargets->find(targJ) != excludeTargets->end()) ||
					(excludeTargetChromosomes != NULL && excludeTargetChromosomes->find(targJ.getChr()) != excludeTargetChromosomes->end()) ||
					targLenFails) {
				cerr << "Excluded target " << targJ;
				if (targLenFails)
					cerr << " of length " << targLen;
				cerr << endl;

				excludeTargetIndices->insert(j);
				excludedTargets->insert(targJ);
			}
		}
	}

	if (!excludeSampleIndices->empty() || !excludeTargetIndices->empty()) {
		HMM_PP::DoubleMat* newRdMat = rdMat->deleteRowsAndColumns(excludeSampleIndices, excludeTargetIndices);
		delete rdMat;
		rdMat = newRdMat;
	}
	delete excludeSampleIndices;
	delete excludeTargetIndices;

	return LoadedReadDepths(rdMat, excludedTargets, excludedSamples);
}

XHMM::xhmmInputManager::IntervalSet* XHMM::xhmmInputManager::readIntervalsFromFile(string intervalsFile, const set<char>& excludeLinesStartingWith) {
	HMM_PP::istreamLineReader* intervalsStream = HMM_PP::utils::getIstreamLineReaderFromFile(intervalsFile);
	if (intervalsStream == NULL)
		throw new Exception("Unable to read table from file '" + intervalsFile + "'");

	IntervalSet* intervSet = new IntervalSet();
	while (!intervalsStream->eof()) {
		string* line = new string();
		*intervalsStream >> *line;
		if (line->empty()) {
			delete line;
			continue;
		}
		char firstChar = (*line)[0];
		if (firstChar != NO_EXCLUDE_CHAR && excludeLinesStartingWith.find(firstChar) != excludeLinesStartingWith.end()) {
			delete line;
			continue;
		}
		stringstream* lineStream = new stringstream(*line);
		delete line;

		intervSet->insert(Interval(*lineStream));
		delete lineStream;
	}
	delete intervalsStream;

	return intervSet;
}

list<string>* XHMM::xhmmInputManager::readStringsFromFile(string stringsFile, char excludeLinesStartingWith) {
	HMM_PP::istreamLineReader* stringsStream = HMM_PP::utils::getIstreamLineReaderFromFile(stringsFile);
	if (stringsStream == NULL)
		throw new Exception("Unable to read from file '" + stringsFile + "'");

	list<string>* stringList = new list<string>();
	while (!stringsStream->eof()) {
		string* line = new string();
		*stringsStream >> *line;
		if (line->empty()) {
			delete line;
			continue;
		}
		char firstChar = (*line)[0];
		if (firstChar != NO_EXCLUDE_CHAR && firstChar == excludeLinesStartingWith) {
			delete line;
			continue;
		}
		stringstream* lineStream = new stringstream(*line);
		delete line;

		while (*lineStream && !lineStream->eof()) {
			string str;
			// Instead of splitting by whitespace, use ONLY tab as delimiter for the header line (to allow for sample names with space in them):
			if (!getline(*lineStream, str, '\t'))
				break;

			stringList->push_back(str);
		}
		delete lineStream;
	}
	delete stringsStream;

	return stringList;
}

XHMM::xhmmInputManager::OneToOneStrMap* XHMM::xhmmInputManager::tableToMap(Table* table, ullint fromColumn, ullint toColumn) {
	if (table == NULL)
		return NULL;

	const ullint ncol = table->getNumColumns();
	if (fromColumn > ncol || toColumn > ncol || fromColumn == 0 || toColumn == 0) {
		stringstream str;
		str << "Cannot extract columns " << fromColumn << "," << toColumn << " from table.";
		throw new Exception(str.str());
	}

	// Access to 0-based indices:
	--fromColumn;
	--toColumn;

	OneToOneStrMap* oom = new OneToOneStrMap();
	for (ullint i = 0; i < table->getNumRows(); ++i)
		oom->addMapping(table->getEntry(i, fromColumn), table->getEntry(i, toColumn));
	delete table;

	return oom;
}

set<string>* XHMM::xhmmInputManager::tableColumnToSet(Table* table, ullint column) {
	if (table == NULL)
		return NULL;

	const ullint ncol = table->getNumColumns();
	if (column > ncol || column == 0) {
		stringstream str;
		str << "Cannot extract column " << column << " from table.";
		throw new Exception(str.str());
	}

	// Access to 0-based indices:
	--column;

	set<string>* vals = new set<string>();
	for (ullint i = 0; i < table->getNumRows(); ++i)
		vals->insert(table->getEntry(i, column));

	delete table;

	return vals;
}
