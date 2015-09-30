#ifndef __DATA_LOADER_H__
#define __DATA_LOADER_H__

#include "utils.hpp"
#include "Data.hpp"
#include "Exception.hpp"
#include "utils.hpp"

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

namespace HMM_PP {

	template<class DataType>
	class DataLoader {

	public:
		DataLoader(string dataFile, bool readHeaderLine = true, bool readDataIDs = true, const set<string>* keepIDs = NULL);
		virtual ~DataLoader();

		virtual bool hasNext();
		virtual DataType* next();

		const string* getHeader() const { return _header; }

	protected:
		istreamLineReader* _stream;

		string* _header;

		bool _readDataIDs;
		const set<string>* _keepIDs;
		uint _idCount;
		uint _numRead;

		DataType* _nextToRead;

		virtual void setValFromStream(typename DataType::InputType& f, istream& stream);

		virtual bool hasMoreStreaming();
		virtual void advanceToNext();
	};

	template<class DataType>
	DataLoader<DataType>::DataLoader(string dataFile, bool readHeaderLine, bool readDataIDs, const set<string>* keepIDs)
	: _stream(utils::getIstreamLineReaderFromFile(dataFile)), _header(NULL), _readDataIDs(readDataIDs), _keepIDs(keepIDs), _idCount(0), _numRead(0), _nextToRead(NULL) {

		if (readHeaderLine) {
			if (!hasMoreStreaming())
				throw new Exception("Unable to read header line from input stream");

			_header = new string();
			*_stream >> *_header;

			if (readDataIDs) {
				stringstream* lineStream = new stringstream(*_header);

				string dummy;
				// Instead of splitting by whitespace, use ONLY tab as delimiter for the row ID ("sample name"):
				if (!getline(*lineStream, dummy, '\t') || !*lineStream)
					throw new Exception("Data input stream failed while reading Matrix title");

				// Read the remainder of the line into _header:
				stringbuf* headerBuf = new stringbuf();
				lineStream->get(*headerBuf);
				if (!*lineStream)
					throw new Exception("Unable to read header line from input stream");
				delete lineStream;
				*_header = headerBuf->str();
				delete headerBuf;
			}
		}
	}

	template<class DataType>
	DataLoader<DataType>::~DataLoader() {
		if (_stream != NULL)
			delete _stream;

		if (_header != NULL)
			delete _header;
	}
	
	template<class DataType>
	bool DataLoader<DataType>::hasMoreStreaming() {
		return !_stream->eof();
	}

	template<class DataType>
	bool DataLoader<DataType>::hasNext() {
		if (_nextToRead == NULL)
			advanceToNext();

		return (_nextToRead != NULL); // still something to be returned in next()
	}

	template<class DataType>
	void DataLoader<DataType>::advanceToNext() {
		// Reset _nextToRead to NULL:
		_nextToRead = NULL;

		if (_keepIDs != NULL && _numRead == _keepIDs->size())
			// In practice, can never read more rows than those in _keepIDs. So, if already found all of _keepIDs, can exit early.
			return;

		while (_nextToRead == NULL) {
			if (!hasMoreStreaming())
				return;

			string* line = new string();
			*_stream >> *line;
			stringstream* lineStream = new stringstream(*line);
			delete line;

			string id;
			if (_readDataIDs) {
				// Instead of splitting by whitespace, use ONLY tab as delimiter for the row ID ("sample name"):
				if (!getline(*lineStream, id, '\t') || !*lineStream)
					throw new Exception("Unable to read data ID from stream");
			}
			else {
				stringstream str;
				str << "ID" << ++_idCount;
				id = str.str();
			}

			// Either no list given, or id is on the list to keep:
			if (_keepIDs == NULL || _keepIDs->find(id) != _keepIDs->end()) {
				++_numRead;

				_nextToRead = new DataType();
				_nextToRead->setId(id);

				while (*lineStream && !lineStream->eof()) {
					typename DataType::InputType f;
					setValFromStream(f, *lineStream);
					_nextToRead->addDatapoint(f);
				}
			}

			delete lineStream;
		}
	}

	template<class DataType>
	DataType* DataLoader<DataType>::next() {
		if (!hasNext())
			throw new Exception("Illegal call to next() when !hasNext()");

		DataType* d = _nextToRead;
		_nextToRead = NULL;

		return d;
	}

	template<class DataType>
	void DataLoader<DataType>::setValFromStream(typename DataType::InputType& f, istream& stream) {
		stream >> f;
		if (!stream)
			throw new Exception("Data input stream failed while reading value");
	}
}

#endif
