#ifndef __ON_DISK_MATRIX_TRANSPOSER_H__
#define __ON_DISK_MATRIX_TRANSPOSER_H__

#include <params.hpp>

#include <ostream>
#include <vector>
#include <set>
#include <map>
#include <list>
using namespace std;

namespace XHMM {
	class OnDiskMatrixTransposer {

	public:
		OnDiskMatrixTransposer(string outFile, const uint numCols);
		virtual ~OnDiskMatrixTransposer();

		void cacheNextRowToDBandDelete(const list<string>* rowData);

		class PullOutputData {
		public:
			PullOutputData() {}
			virtual ~PullOutputData() {}

			virtual void setNextDataPoint(const string& data, int ind) = 0;
		};
		void nextTransposedRowFromDB(PullOutputData* pullData);

	private:
		void openDBfileWriters();
		void closeDBfileWriters();

		void deleteDBfiles();

		void openDBfileReaders();
		void closeDBfileReaders();

		uint _numRows;
		const uint _numCols;

		string _dbFile;

		ostream* _dbFileWriter;
		ostream* _dbIndexFileWriter;

		list<istream*>* _columnReaders;
	};

}

#endif
