#include "utils.hpp"
#include "Exception.hpp"

#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>
using namespace std;

//////////////////////////////////////////////////////////////
//Implementation of HMM_PP::istreamLineReader
//////////////////////////////////////////////////////////////
const string HMM_PP::istreamLineReader::BLANK = "";

HMM_PP::istreamLineReader::istreamLineReader(istream& stream, bool deleteInner)
: _innerStream(stream), _deleteInner(deleteInner), _nextToExtract(BLANK) {
}

HMM_PP::istreamLineReader::~istreamLineReader() {
	if (_deleteInner)
		delete &_innerStream;
}

bool HMM_PP::istreamLineReader::eof() {
	if (_nextToExtract == BLANK)
		advanceToNonBlankLine();

	return (_nextToExtract == BLANK); //nothing left to extract
}

HMM_PP::istreamLineReader& HMM_PP::istreamLineReader::operator>>(string& s) {
	if (this->eof()) // since the user is requesting to read input, this advances to the next line ["just-in-time" paradigm]
		return *this;

	s = _nextToExtract;
	_nextToExtract = BLANK;

	return *this;
}

void HMM_PP::istreamLineReader::advanceToNonBlankLine() {
	/* Reset _nextToExtract to a non-extractable value of BLANK
     [since if _nextToExtract were truly BLANK,
     then only have a '\n', but then will have skipped such a line]: */
	_nextToExtract = BLANK;
	if (!_innerStream.good()) //shouldn't do peek() if !good()
		return;

	//Skip '\n' lines:
	for (char nextChar = _innerStream.peek(); !_innerStream.eof(); nextChar = _innerStream.peek()) {
		if (nextChar == '\n') //blank line, so continue [otherwise get(lineBuf) would fail]
			_innerStream.get(); //get '\n'
		else
			break;
	}
	if (_innerStream.eof()) //shouldn't do get() if eof()
		return;

	stringbuf lineBuf;
	_innerStream.get(lineBuf);
	_innerStream.get(); //get '\n'
	_nextToExtract = lineBuf.str();

	stringstream nextLineStream;
	nextLineStream << _nextToExtract;
	bool nonWhiteSpace = false;
	while (!nextLineStream.eof()) {
		string tmpStr;
		nextLineStream >> tmpStr;
		if (tmpStr != BLANK) {
			nonWhiteSpace = true;
			break;
		}
	}
	if (!nonWhiteSpace) //recursively get rid of '\n' lines and/or white-space-only lines
		advanceToNonBlankLine();
}
//////////////////////////////////////////////////////////////
//End Implementation of HMM_PP::istreamLineReader
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//Implementation of HMM_PP::ostreamWriter
//////////////////////////////////////////////////////////////
HMM_PP::ostreamWriter::ostreamWriter(ostream& stream, bool deleteInner)
: _stream(stream), _deleteInner(deleteInner) {
}

HMM_PP::ostreamWriter::~ostreamWriter() {
	if (_deleteInner)
		delete &_stream;
}
//////////////////////////////////////////////////////////////
//End Implementation of HMM_PP::ostreamWriter
//////////////////////////////////////////////////////////////

const string HMM_PP::utils::STD_STREAM = "-";

HMM_PP::istreamLineReader* HMM_PP::utils::getIstreamLineReaderFromFile(string file, bool requireNamedFile) {
	if (file == "")
		return NULL;

	if (file == STD_STREAM) {
		if (requireNamedFile)
			throw new Exception("Cannot open standard input as a named file");
		return new istreamLineReader(cin);
	}

	istream* stream = new ifstream(file.c_str());
	if (!*stream || stream->eof())
		throw new FileNotFoundException("File " + file + " does not exist, could not be opened for reading, or is empty");
	return new istreamLineReader(*stream, true);
}

HMM_PP::ostreamWriter* HMM_PP::utils::getOstreamWriterFromFile(string file) {
	if (file == "")
		return NULL;

	if (file == STD_STREAM)
		return new ostreamWriter(cout);

	ofstream* stream = new ofstream(file.c_str());
	if (!*stream)
		throw new Exception("File " + file + " could not be opened for writing");
	return new ostreamWriter(*stream, true);
}


// Helper functions
string HMM_PP::int2str(const uint i)  {
	stringstream ss; ss << i;
	return ss.str();
}

bool HMM_PP::fileExists(const string& f) {
	return access(f.c_str(), F_OK) == 0;
}

bool HMM_PP::removeFile(const string& f) {
	return remove(f.c_str()) == 0;
}

bool HMM_PP::beginsWith(string str, string prefix) {
	if (str.size() < prefix.size())
		return false;

	return equal(prefix.begin(), prefix.end(), str.begin());
}

bool HMM_PP::endsWith(string str, string suffix) {
	if (str.size() < suffix.size())
		return false;

	return equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

// Helper functions
bool HMM_PP::realNum(double d) {
	double zero = 0;
	if (d != d || d == 1/zero || d == -1/zero)
		return false;
	else
		return true;
}

/*
double HMM_PP::sumLogProb(double logprob1, double logprob2) {
	if (isinf(logprob1)&& isinf(logprob2))
		return logprob1; // both prob1 and prob2 are 0, return log 0.
	if (logprob1>logprob2)
		return logprob1+log(1+exp(logprob2-logprob1));
	else
		return logprob2+log(1+exp(logprob1-logprob2));
}

double HMM_PP::sumLogProb(const vector<double>& logprobs) {
	double max = 0;
	uint i;
	for (i = 0; i<logprobs.size(); i++) {
		if (i==0 || logprobs[i]>max)
			max = logprobs[i];
	}
	if (isinf(max)) // the largest probability is 0 (log prob= -inf)
		return max;   // return log 0
	double p = 0;
	for (i = 0; i<logprobs.size(); i++) {
		p += exp(logprobs[i]-max);
	}
	return max + log(p);
}
 */

/*
bool HMM_PP::constrain(double* v, double l, double u) {
	bool inrange = true;
	if (*v < l)  {
		inrange = false;
 *v = l;
	}
	else if (*v > u) {
		inrange = false;
 *v = u;
	}
	return inrange;
}
 */
