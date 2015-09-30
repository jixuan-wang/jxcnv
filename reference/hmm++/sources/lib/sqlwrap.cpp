#include "sqlwrap.hpp"

bool SQL::open(const string& n, const string& scratch) {

	_name = n;

	_rc = sqlite3_open(_name.c_str(),&_db);

	if (_rc) halt("problem opening database: " + _name);

	if (scratch  != "") {
		query("PRAGMA temp_store_directory = '"
				+ scratch
				+ "';");
	}
	return _rc == 0;
}

void SQL::synchronous(bool b) {
	if (!b)
		query("PRAGMA synchronous=OFF;");
	else
		query("PRAGMA synchronous=FULL;");
}

bool SQL::table_exists(const string& table_name) {
	sqlite3_stmt* s = prepare("SELECT name FROM sqlite_master WHERE type='table' AND name= :table_name ; ");
	bind_text(s, ":table_name", table_name);
	if (step(s))  {
		finalize(s);
		return true;
	}
	finalize(s);
	return false;
}

void SQL::close() {
	if (_db)  {
		sqlite3_close(_db);
		_db = NULL;
	}
}

bool SQL::query(const string& q) {
	char* db_err;
	_rc = sqlite3_exec(_db, q.c_str(), 0, 0,&db_err);
	if (_rc) warn(string(db_err));
	return _rc == 0;
}

sqlite3_stmt* SQL::prepare(const string& q) {
	sqlite3_stmt* p;
	int rc = sqlite3_prepare_v2(_db, q.c_str(), q.size(),&p, NULL);
	if (rc) warn("preparing query " + string(sqlite3_errmsg(_db)));
	else _qset.insert(p);
	return rc ? NULL : p;
}

sqlite3_stmt* SQL::prepare(const string& q, const string& key) {
	sqlite3_stmt* p;
	int rc = sqlite3_prepare(_db, q.c_str(), q.size(),&p, NULL);
	if (rc) halt(_db_err);
	else _qset.insert(p);
	_qmap.insert(make_pair(key, p));
	return rc ? NULL : p;
}

sqlite3_stmt* SQL::fetch_prepared(const string& key) {
	map<string,sqlite3_stmt*>::iterator i = _qmap.find(key);
	if (i == _qmap.end()) return NULL;
	return i->second;
}

void SQL::begin() {
	char* db_err;
	string q = "BEGIN;";
	_rc = sqlite3_exec(_db, q.c_str(), 0, 0,&db_err);
	if (_rc) halt(db_err);
}

void SQL::finalize(sqlite3_stmt* stmt) {
	set<sqlite3_stmt*>::iterator i = _qset.find(stmt);
	if (stmt&& i != _qset.end())  {
		_qset.erase(i);
		sqlite3_finalize(stmt);
	}
	stmt = NULL;
}

bool SQL::step(sqlite3_stmt* stmt) {
	_rc = sqlite3_step(stmt);

	if (_rc != SQLITE_ROW && _rc != SQLITE_DONE) {
		reset(stmt);
		halt(sqlite3_errmsg(_db));
	}

	return _rc == SQLITE_ROW;
}

void SQL::reset(sqlite3_stmt* stmt) {
	sqlite3_reset(stmt);
}


void SQL::bind_int(sqlite3_stmt* stmt, const string& index, int value) {
	sqlite3_bind_int(stmt,
			sqlite3_bind_parameter_index(stmt, index.c_str()),
			value);
}

void SQL::bind_null(sqlite3_stmt* stmt, const string& index ) {
	sqlite3_bind_null(stmt,
			sqlite3_bind_parameter_index(stmt, index.c_str()));
}

void SQL::bind_int64(sqlite3_stmt* stmt, const string& index, sqlite_int64 value) {
	sqlite3_bind_int64(stmt,
			sqlite3_bind_parameter_index(stmt, index.c_str()),
			value);
}

void SQL::bind_double(sqlite3_stmt* stmt, const string& index, double value) {
	sqlite3_bind_double(stmt,
			sqlite3_bind_parameter_index(stmt, index.c_str()),
			value);
}

void SQL::bind_text(sqlite3_stmt* stmt, const string& index, const string& value) {

	sqlite3_bind_text(stmt,
			sqlite3_bind_parameter_index(stmt, index.c_str()),
			value.c_str(),
			value.size(),
			0);
}


void SQL::bind_blob(sqlite3_stmt* stmt, const string& index, blob& value) {
	_rc = sqlite3_bind_blob(stmt,
			sqlite3_bind_parameter_index(stmt, index.c_str()),
			value._p,
			value._l,
			0);

}


int SQL::get_int(sqlite3_stmt* stmt, int idx) {
	return sqlite3_column_int(stmt, idx);
}

sqlite_int64 SQL::get_int64(sqlite3_stmt* stmt, int idx) {
	return sqlite3_column_int64(stmt, idx);
}

double SQL::get_double(sqlite3_stmt* stmt, int idx) {
	return sqlite3_column_double(stmt, idx);
}

string SQL::get_text( sqlite3_stmt* stmt, int idx) {
	const unsigned char* s = sqlite3_column_text(stmt, idx);
	if (s == NULL)
		return "";
	else
		return (const char*)s;
}

blob SQL::get_blob(sqlite3_stmt* stmt, int idx) {
	blob b;
	b._p = (const char*)sqlite3_column_blob(stmt, idx);
	b._l = sqlite3_column_bytes(stmt, idx);
	return b;
}

void SQL::commit() {
	query("COMMIT;");
}

vector<int> SQL::intTable(const string& q, int cols) {
	return intTable(prepare(q), cols);
}

vector<int> SQL::intTable(sqlite3_stmt* stmt, int cols) {
	vector<int> res;
	_rc = sqlite3_step(stmt);
	while (_rc == SQLITE_ROW) {
		for (int i = 0 ; i < cols ; i++)
			res.push_back (sqlite3_column_int(stmt, i));
		_rc = sqlite3_step(stmt);
	}
	sqlite3_finalize(stmt);
	return res;
}

vector<sqlite_int64> SQL::int64Table(const string& q, int cols) {
	return int64Table(prepare(q), cols);
}

vector<sqlite_int64> SQL::int64Table(sqlite3_stmt* stmt, int cols) {
	vector<sqlite_int64> res;
	_rc = sqlite3_step(stmt);
	while (_rc == SQLITE_ROW) {
		for (int i = 0 ; i < cols ; i++)
			res.push_back (sqlite3_column_int64(stmt, i));
		_rc = sqlite3_step(stmt);
	}
	sqlite3_finalize(stmt);
	return res;
}

int SQL::lookup_int(sqlite3_stmt* stmt) {
	int r = -1;
	_rc = sqlite3_step(stmt);
	if (_rc == SQLITE_ROW)
		r = sqlite3_column_int(stmt, 0);
	return r;
}

int SQL::lookup_int(const string& q) {
	sqlite3_stmt* s = prepare(q);
	int r = -1;
	_rc = sqlite3_step(s);
	if (_rc == SQLITE_ROW)
		r = sqlite3_column_int(s, 0);
	finalize(s);
	return r;
}

sqlite_int64 SQL::lookup_int64(sqlite3_stmt* stmt) {
	sqlite_int64 r = 0;
	_rc = sqlite3_step(stmt);
	if (_rc == SQLITE_ROW)
		r = sqlite3_column_int64(stmt, 0);
	return r;
}

