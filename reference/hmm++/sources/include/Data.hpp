#ifndef __DATA_H__
#define __DATA_H__

#include "params.hpp"

#include <string>
#include <vector>
#include <iostream>
using namespace std;

class SQL;
struct sqlite3_stmt;

namespace HMM_PP {
	// Observations for a particular sequence
	class Data {

	public:
		Data();
		virtual ~Data() = 0;

	protected:
		virtual void clearObservations() = 0;

	public:
		void clear();

		const string& getId() const { return _id; }
		void setId(const string& i) { _id = i; }

		virtual inline void setNumObservations(const uint nobs) { _n = nobs; }
		virtual inline uint getNumObservations() const { return _n; }

		virtual void printDatapoint(ostream& stream, const uint t) const = 0;

		virtual void bindDatapointToDB(SQL& db, sqlite3_stmt* stmt, const string& data, const uint t) = 0;
		virtual void addDatapointFromDB(SQL& db, sqlite3_stmt* stmt, const uint queryPos) = 0;

		uint checkStartStopRange(const uint t1, const vector<uint>& statePath) const;
		void checkStartStopRange(const uint t1, const uint t2) const;

	protected:
		uint _n;
		string _id;
	};
}

#endif
