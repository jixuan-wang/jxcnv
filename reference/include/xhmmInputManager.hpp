#ifndef __XHMM_INPUT_MANAGER_H__
#define __XHMM_INPUT_MANAGER_H__

#include "Interval.hpp"

#include <Exception.hpp>
#include <params.hpp>
#include <utils.hpp>
#include <UnivariateDataNoDB.hpp>
#include <DataLoader.hpp>

#include <vector>
#include <string>
#include <map>
#include <list>
#include <set>
#include <sstream>
using namespace std;

#define NO_EXCLUDE_CHAR ' '

#define COMMENT_LINE_CHAR '#'
#define INTERVALS_HEADER_LINE_CHAR '@'

namespace HMM_PP {
	template<class T> class NamedMatrix;
}

namespace XHMM {

	class xhmmInputManager {

	private:
		static set<char> createIntervalsExcludeCharSet();

	public:
		class Table {
		public:
			friend class xhmmInputManager;
			typedef vector<string> Row;

			Table(Row* colNames);

		private:
			Table();

		public:
			~Table();

			const Row& getColumns() const { return *_colNames; }
			ullint getNumColumns() const { return _colNames->size(); }
			bool hasColumn(string name) const { return _colNameToIndex->find(name) != _colNameToIndex->end(); }

			ullint getNumRows() const { return _rows->size(); }

			const string& getEntry(ullint row, const string& column) const;
			const string& getEntry(ullint row, ullint column) const;

			void addRow(Row* row);

		private:
			Row* _colNames;
			map<string, ullint>* _colNameToIndex;

			vector<Row*>* _rows;
		};

		static Table* readTable(string tableFile, bool header = true);
		static Table* readFastaIndexTable(string fastaFile);

		static list<string>* getListFromArray(char** vals, const ullint length);

		typedef set<Interval> IntervalSet;
		typedef set<string> StringSet;

		template<class BaseType>
		class MatrixReader {
		public:
			static HMM_PP::NamedMatrix<BaseType>* readMatrixFromFile(string file, bool colnames = true, bool rownames = true, ullint maxNumRows = 0);
		};

		class LoadedReadDepths {
		public:
			LoadedReadDepths(HMM_PP::DoubleMat* rd, IntervalSet* excludedTargets, StringSet* excludedSamples)
			: _rd(rd), _excludedTargets(excludedTargets), _excludedSamples(excludedSamples) {}
			~LoadedReadDepths() {}

			HMM_PP::DoubleMat* _rd;
			IntervalSet* _excludedTargets;
			StringSet* _excludedSamples;
		};

		static LoadedReadDepths loadReadDepthsFromFile(string readDepthFile, IntervalSet* excludeTargets = NULL, StringSet* excludeTargetChromosomes = NULL, StringSet* excludeSamples = NULL,
				const ullint minTargetSize = 0, const ullint maxTargetSize = ULLINT_INFINITY);

		static set<char> INTERVALS_EXCLUDE_LINE_CHARS;
		static IntervalSet* readIntervalsFromFile(string intervalsFile, const set<char>& excludeLinesStartingWith = INTERVALS_EXCLUDE_LINE_CHARS);

		static list<string>* readStringsFromFile(string stringsFile, char excludeLinesStartingWith = COMMENT_LINE_CHAR);

		template<class From, class To>
		class OneToOneMap {
		public:
			OneToOneMap() {}
			~OneToOneMap() {}

			void addMapping(From f, To t) {
				if (_fromToMap.find(f) != _fromToMap.end())
					throw new Exception("Cannot add mapping from " + f + " more than once.");
				if (_toSet.find(t) != _toSet.end())
					throw new Exception("Cannot add mapping to " + t + " more than once.");

				_fromToMap[f] = t;
				_toSet.insert(t);
			}

			typedef map<From, To> FromToMap;
			typedef typename FromToMap::const_iterator iter;

			iter find(const From& f) { return _fromToMap.find(f); }
			iter end() { return _fromToMap.end(); }

			const To& operator[](const From& f) const { return _fromToMap[f]; }

		private:
			FromToMap _fromToMap;
			set<To> _toSet;
		};

		typedef OneToOneMap<string, string> OneToOneStrMap;
		static OneToOneStrMap* tableToMap(Table* table, ullint fromColumn, ullint toColumn);

		static set<string>* tableColumnToSet(Table* table, ullint column = 1);

	private:
		static Table::Row* readRowFromStream(HMM_PP::istreamLineReader* tableStream);

		static const string FASTA_SUFFIX;
		static const string FASTA_INDEX_SUFFIX;
	};

	template<class BaseType>
	HMM_PP::NamedMatrix<BaseType>* XHMM::xhmmInputManager::MatrixReader<BaseType>::readMatrixFromFile(string file, bool colnames, bool rownames, ullint maxNumRows) {
		typedef HMM_PP::UnivariateDataNoDB<BaseType> DataType;

		HMM_PP::DataLoader<DataType>* loader = new HMM_PP::DataLoader<DataType>(file, colnames, rownames);

		ullint numColumns = 0;
		list<string>* columns = NULL;
		if (colnames) {
			columns = new list<string>();
			stringstream* lineStream = new stringstream(*(loader->getHeader()));

			while (*lineStream && !lineStream->eof()) {
				string colString;
				// Instead of splitting by whitespace, use ONLY tab as delimiter for the column headers:
				if (!getline(*lineStream, colString, '\t') || !*lineStream)
					throw new Exception("Data input stream failed while reading column " + colString);
				columns->push_back(colString);
			}
			delete lineStream;

			numColumns = columns->size();
		}

		HMM_PP::NamedMatrix<BaseType>* mat = new HMM_PP::NamedMatrix<BaseType>(0, columns->size());

		ullint row = 0;
		while (loader->hasNext() && (maxNumRows == 0 || row < maxNumRows)) {
			DataType* data = loader->next();
			const ullint numObs = data->getNumObservations();

			if (numColumns == 0)
				numColumns = numObs;
			else if (numObs != numColumns) {
				stringstream str;
				str << "Cannot read jagged matrix of unequal number of columns - expected: " << numColumns << " found: " << numObs;
				throw new Exception(str.str());
			}

			mat->addRow(data->getValuesTransferOwner());
			if (rownames)
				mat->setRowName(row, data->getId());

			delete data;
			++row;
		}

		delete loader;

		if (colnames) {
			ullint col = 0;
			for (list<string>::const_iterator colIt = columns->begin(); colIt != columns->end(); ++colIt)
				mat->setColName(col++, *colIt);
			delete columns;
		}

		return mat;
	}

}

#endif
