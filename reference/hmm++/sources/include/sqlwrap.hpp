#ifndef __SQLWRAP_H__
#define __SQLWRAP_H__

#include "sqlite3.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cstdlib>
using namespace std;

class blob {

public:
	string _s;
	const char* _p;
	int _l;

	blob() {}

	blob(const string& t)  {
		_s = t;
		_p = _s.data();
		_l = _s.size();
	}

	void set_string(const string& tmp) {
		_s = tmp;
		_p = _s.data();
		_l = _s.size();
	}

	string get_string() {
		return string((const char*)_p,_l);
	}
};


class SQL {

public:

	SQL() {
		_db = NULL;
		_name = "";
	}

	bool open(const string& n, const string& scratch = "");

	void synchronous(bool);
	void close();
	bool attached() const { return _db; }
	bool is_open() const { return _db; }
	bool query(const string& q);
	bool table_exists(const string&);

	void halt(const string& msg)  {
		cerr << "_db error : " << msg << endl;
		exit(1);
	}

	void warn(const string& msg)  {
		cerr << "_db warning : " << msg << endl;
	}

	sqlite3_stmt* prepare(const string& q);
	sqlite3_stmt* prepare(const string& q, const string& key);
	sqlite3_stmt* fetch_prepared(const string& key);

	bool step(sqlite3_stmt* stmt);
	void reset(sqlite3_stmt* stmt);
	void finalize(sqlite3_stmt* stmt);

	bool loadExtension(string);

	void begin();
	void commit();

	sqlite_int64 last_insert_rowid() { return sqlite3_last_insert_rowid(_db); }

	void bind_int(sqlite3_stmt* stmt, const string& index, int value);
	void bind_int64(sqlite3_stmt* stmt, const string& index, sqlite_int64 value);
	void bind_double(sqlite3_stmt* stmt, const string& index, double value);
	void bind_text(sqlite3_stmt* stmt, const string& index, const string& value);
	void bind_blob(sqlite3_stmt* stmt, const string& index, blob&);
	void bind_null(sqlite3_stmt* stmt, const string& index);

	int get_int(sqlite3_stmt*, int);
	sqlite_int64 get_int64(sqlite3_stmt*, int);
	double get_double(sqlite3_stmt*, int);
	string get_text(sqlite3_stmt*, int);
	blob get_blob(sqlite3_stmt*, int);

	int lookup_int(sqlite3_stmt*);
	int lookup_int(const string& q);
	sqlite_int64 lookup_int64(sqlite3_stmt*);
	vector<int> intTable(const string& q, int cols);
	vector<int> intTable(sqlite3_stmt* stmt, int cols);

	vector<sqlite_int64> int64Table(const string& q, int cols);
	vector<sqlite_int64> int64Table(sqlite3_stmt* stmt, int cols);

	static string header_version() {
		return sqlite3_libversion();
	}

	static string library_version() {
		return SQLITE_VERSION;
	}

	sqlite3* pointer() { return _db; }

private:
	// Keep track of all prepared statements
	set<sqlite3_stmt*> _qset;

	// Map of prepared statements
	map<string,sqlite3_stmt*> _qmap;

	// Database
	sqlite3* _db;

	// Return code, error msg
	int _rc;
	char* _db_err;

	// Name of database
	string _name;
};

#endif
