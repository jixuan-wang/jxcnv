#include "OnDiskMatrixTransposer.hpp"
#include "xhmmOutputManager.hpp"

#include <params.hpp>
#include <utils.hpp>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
using namespace std;

#define DB_SUFFIX ".db"
#define DB_INDEX_SUFFIX ".dbi"

XHMM::OnDiskMatrixTransposer::OnDiskMatrixTransposer(string outFile, const uint numCols)
: _numRows(0), _numCols(numCols),
  _dbFile(outFile + DB_SUFFIX),
  _dbFileWriter(NULL), _dbIndexFileWriter(NULL), _columnReaders(NULL) {

	openDBfileWriters();
}

XHMM::OnDiskMatrixTransposer::~OnDiskMatrixTransposer() {
	closeDBfileWriters();
	closeDBfileReaders();
	deleteDBfiles();
}

/*
 * DB functions:
 */
void XHMM::OnDiskMatrixTransposer::openDBfileWriters() {
	closeDBfileWriters();
	deleteDBfiles();

	// open in binary to permit proper tellp() [i.e., handling of new lines]:
	_dbFileWriter = new ofstream(_dbFile.c_str(), ofstream::out | ofstream::trunc | ofstream::binary);
	_dbIndexFileWriter = new ofstream((_dbFile + DB_INDEX_SUFFIX).c_str());

	if (!*_dbFileWriter || !*_dbIndexFileWriter)
		throw new Exception("Unable to open DB file " + _dbFile + " and its index for writing");
}

void XHMM::OnDiskMatrixTransposer::closeDBfileWriters() {
	if (_dbFileWriter != NULL) {
		delete _dbFileWriter;
		_dbFileWriter = NULL;
	}

	if (_dbIndexFileWriter != NULL) {
		delete _dbIndexFileWriter;
		_dbIndexFileWriter = NULL;
	}
}

void XHMM::OnDiskMatrixTransposer::deleteDBfiles() {
	if (HMM_PP::fileExists(_dbFile) && !HMM_PP::removeFile(_dbFile))
		xhmmOutputManager::printWarning("Failed to remove DB file " + _dbFile);

	if (HMM_PP::fileExists(_dbFile + DB_INDEX_SUFFIX) && !HMM_PP::removeFile(_dbFile + DB_INDEX_SUFFIX))
		xhmmOutputManager::printWarning("Failed to remove DB file " + _dbFile + DB_INDEX_SUFFIX);
}

void XHMM::OnDiskMatrixTransposer::openDBfileReaders() {
	closeDBfileWriters();
	closeDBfileReaders();

	ifstream indexStream((_dbFile + DB_INDEX_SUFFIX).c_str());
	if (!indexStream)
		throw new Exception("Unable to open DB index " + _dbFile + DB_INDEX_SUFFIX + " for reading");

	_columnReaders = new list<istream*>();
	for (uint i = 0; i < _numRows; ++i) {
		// open in binary to permit proper seekg() [i.e., handling of new lines]:
		istream* reader = new ifstream(_dbFile.c_str(), ifstream::in | ifstream::binary);
		istream::streamoff nextRowPos;
		indexStream >> nextRowPos;
		if (!indexStream)
			throw new Exception("Unable to read from DB index " + _dbFile + DB_INDEX_SUFFIX);

		reader->seekg(nextRowPos);
		_columnReaders->push_back(reader);
	}
}

void XHMM::OnDiskMatrixTransposer::closeDBfileReaders() {
	if (_columnReaders != NULL) {
		for (list<istream*>::const_iterator it = _columnReaders->begin(); it != _columnReaders->end(); ++it)
			delete *it;
		delete _columnReaders;
		_columnReaders = NULL;
	}
}

void XHMM::OnDiskMatrixTransposer::cacheNextRowToDBandDelete(const list<string>* rowData) {
	if (rowData->size() != _numCols) {
		stringstream str;
		str << "Must provide rows of size " << _numCols << " [!= " << rowData->size() << "]";
		throw new Exception(str.str());
	}

	ostream::streamoff nextRowPos = _dbFileWriter->tellp();
	if (nextRowPos == -1)
		throw new Exception("Unable to get DB write pointer for row");
	*_dbIndexFileWriter << nextRowPos << endl;
	if (!*_dbIndexFileWriter)
		throw new Exception("Unable to write to DB index");

	uint colInd = 0;
	for (list<string>::const_iterator it = rowData->begin(); it != rowData->end(); ++it) {
		const string& dataStr = *it;
		if (colInd++ > 0)
			*_dbFileWriter << '\t';
		*_dbFileWriter << dataStr;
	}
	*_dbFileWriter << endl;
	if (!*_dbFileWriter)
		throw new Exception("Unable to write to DB");

	delete rowData;

	++_numRows;
}

void XHMM::OnDiskMatrixTransposer::nextTransposedRowFromDB(PullOutputData* pullData) {
	if (_columnReaders == NULL) {
		cerr << "Preparing to read transposed rows from disk..." << endl;
		openDBfileReaders();
	}

	int ind = 0;
	for (list<istream*>::const_iterator readIt = _columnReaders->begin(); readIt != _columnReaders->end(); ++readIt) {
		istream* reader = *readIt;
		string dataStr;
		*reader >> dataStr;
		if (!*reader) {
			stringstream str;
			str << "Unable to read data from DB row " << ind << " (perhaps too many file handles?)";
			throw new Exception(str.str());
		}

		pullData->setNextDataPoint(dataStr, ind);
		++ind;
	}
}
