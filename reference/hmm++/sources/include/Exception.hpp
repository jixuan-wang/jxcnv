#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

#include <iostream>
#include <string>
using namespace std;

class Exception {
public:

	Exception(string msg) : _msg(msg) {}
	virtual ~Exception() {}

	const string& getMessage() const { return _msg;}
	void printErrorMessage() const { cerr << endl << "ERROR: " << _msg << endl;}

protected:
	string _msg;
};


class FileNotFoundException : public Exception {
public:

	FileNotFoundException(string msg) : Exception(msg) {}
	virtual ~FileNotFoundException() {}
};


#endif
