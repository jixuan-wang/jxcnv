#ifndef __NAMED_MATRIX_H__
#define __NAMED_MATRIX_H__

#include "NamedVector.hpp"
#include "utils.hpp"

#include <vector>
#include <iomanip>
#include <set>
#include <sstream>
using namespace std;

#define CONST_LENGTH_ROWS true

namespace HMM_PP {

	template<class RealType>
	class NamedMatrix {

	public:
		typedef NamedVector<RealType> INNER;
		typedef vector<INNER*> PARENT;

		explicit NamedMatrix(string name)
		: _mat(new PARENT()), _ncol(0), _defaultValue(RealType()), _name(name), _rowNames(NULL), _colNames(NULL) {
		}

		explicit NamedMatrix(const ullint nrow = 0, const ullint ncol = 0, const RealType val = RealType(), string name = "")
		: _mat(new PARENT()), _ncol(0), _defaultValue(RealType()), _name(name), _rowNames(NULL), _colNames(NULL) {
			setDims(nrow, ncol, val);
		}

		NamedMatrix(const NamedVector<RealType>& vec)
		: _mat(new PARENT(1)), _ncol(vec.size()), _defaultValue(RealType()), _name(""), _rowNames(NULL), _colNames(NULL) {
			(*_mat)[0] = new INNER(vec, CONST_LENGTH_ROWS);
		}

	public:
		virtual ~NamedMatrix() {
			this->clearData();
			delete _mat;

			if (_rowNames != NULL)
				delete _rowNames;
			if (_colNames != NULL)
				delete _colNames;
		}

		void addRow(NamedVector<RealType>* row) {
			if (_mat->empty())
				_ncol = row->size();
			else if (row->size() != _ncol)
				throw new Exception("Matrix rows must have uniform number of columns");

			row->setConstantLength();
			_mat->push_back(row);
		}

		void setMatrixName(string name) { _name = name; }
		string getMatrixName() const { if (_name == "") return "Matrix"; return _name; }

		const RealType& operator()(const ullint i, const ullint j) const { return (*this)[i][j]; }
		RealType& operator()(const ullint i, const ullint j) { return (*this)[i][j]; }

		const NamedVector<RealType>& operator[](const ullint i) const { return getRow(i); }
		NamedVector<RealType>& operator[](const ullint i) { return getRow(i); }

		void dropLastRow() {
			if (empty())
				throw new Exception("Cannot delete last row of empty matrix");

			int i = nrow() - 1;
			if (_rowNames != NULL && _rowNames->size() >= nrow())
				_rowNames->resize(i);

			delete (*_mat)[i];
			_mat->resize(i);
		}

	private:
		NamedVector<RealType>& getRow(const ullint i) const { ensureRow(i); return (*(*_mat)[i]); }

		void clearData() {
			for (ullint i = 0; i < _mat->size(); ++i)
				if ((*_mat)[i] != NULL)
				delete (*_mat)[i];

			_mat->clear();
		}

		void turnOnNames() {
			ullintPair dims = getDims();

			if (dims.first > 0) {
				if (_rowNames == NULL)
					_rowNames = new vector<string>(dims.first);
				else
					_rowNames->resize(dims.first);
			}
			else if (_rowNames == NULL)
				_rowNames = new vector<string>();

			if (dims.second > 0) {
				if (_colNames == NULL)
					_colNames = new vector<string>(dims.second);
				else
					_colNames->resize(dims.second);
			}
			else if (_colNames == NULL)
				_colNames = new vector<string>();
		}

		bool namesAreOn() const {
			return _rowNames != NULL && _colNames != NULL;
		}

	public:
#define CHECK_ROW_INDEX \
		if (i > nrow()) { \
			stringstream str; \
			str << "Row index " << i << " is out of bounds"; \
			throw new Exception(str.str()); \
		}

		void setRowName(const ullint i, const string& s) {
			CHECK_ROW_INDEX

			if (!namesAreOn())
				turnOnNames();

			ullint minSz = i + 1;
			if (_rowNames->size() < minSz)
				_rowNames->resize(nrow());

			(*_rowNames)[i] = s;
		}

		string rowName(const ullint i) const {
			CHECK_ROW_INDEX

			ullint minSz = i + 1;
			if (_rowNames == NULL || _rowNames->size() < minSz) {
				stringstream str;
				str << "[" << minSz << ",]";
				return str.str();
			}

			return (*_rowNames)[i];
		}

#define CHECK_COL_INDEX \
		if (j > ncol()) { \
			stringstream str; \
			str << "Col index " << j << " is out of bounds"; \
			throw new Exception(str.str()); \
		}

		void setColName(const ullint j, const string& s) {
			CHECK_COL_INDEX

			if (!namesAreOn())
				turnOnNames();

			ullint minSz = j + 1;
			if (_colNames->size() < minSz)
				_colNames->resize(ncol());

			(*_colNames)[j] = s;
		}

		string colName(const ullint j) const {
			CHECK_COL_INDEX

			ullint minSz = j + 1;
			if (_colNames == NULL || _colNames->size() < minSz) {
				stringstream str;
				str << "[," << minSz << "]";
				return str.str();
			}

			return (*_colNames)[j];
		}

		void setDims(const ullint nrow, const ullint ncol, RealType value = RealType()) {
			this->clearData();

			_ncol = ncol;
			_mat->resize(nrow, NULL);

			_defaultValue = value;
		}

	private:
		// Since _mat is mutable, it can still be lazily set here in this const member function:
		void ensureRow(const ullint i) const {

			if (i >= nrow())
				_mat->resize(i+1, NULL);

			if ((*_mat)[i] == NULL)
				(*_mat)[i] = new INNER(_ncol, _defaultValue, CONST_LENGTH_ROWS);
		}

	public:
		bool empty() const {
			return _mat->empty();
		}

		ullint nrow() const {
			return _mat->size();
		}

		ullint ncol() const {
			return _ncol;
		}

		ullintPair getDims() const {
			return ullintPair(nrow(), ncol());
		}

		static const string TRANSPOSE_SUFFIX;

		NamedMatrix* transpose() const {
			ullintPair rowCol = getDims();
			ullint rows = rowCol.first;
			ullint cols = rowCol.second;

			stringstream nameStr;
			nameStr << getMatrixName() << TRANSPOSE_SUFFIX;
			NamedMatrix* transpose = new NamedMatrix(cols, rows, 0, nameStr.str());
			for (ullint i = 0; i < rows; ++i)
				for (ullint j = 0; j < cols; ++j)
					(*transpose)(j,i) = (*this)(i,j);

			if (_rowNames != NULL)
				transpose->_colNames = new vector<string>(*_rowNames);

			if (_colNames != NULL)
				transpose->_rowNames = new vector<string>(*_colNames);

			return transpose;
		}

		NamedMatrix* deleteRowsAndColumns(const set<ullint>* deleteRows, const set<ullint>* deleteColumns) const {
			ullintPair rowCol = getDims();
			ullint rows = rowCol.first;
			ullint cols = rowCol.second;

			// Filter out any extraneous columns to remove:
			set<ullint>* removeColumns = new set<ullint>();
			if (deleteColumns != NULL) {
				for (set<ullint>::const_iterator colIt = deleteColumns->begin(); colIt != deleteColumns->end(); ++colIt)
					if (*colIt < cols)
						removeColumns->insert(*colIt);
			}

			NamedMatrix* newMat = new NamedMatrix(getMatrixName());
			for (ullint i = 0; i < rows; ++i) {
				if (deleteRows == NULL || deleteRows->find(i) == deleteRows->end()) {
					const INNER* rowI = (*_mat)[i];
					if (removeColumns->empty())
						newMat->addRow(new INNER(*rowI));
					else {
						INNER* newRow = new INNER(rowI->size() - removeColumns->size());
						ullint ind = 0;
						for (ullint j = 0; j < cols; ++j)
							if (removeColumns->find(j) == removeColumns->end())
								(*newRow)[ind++] = (*rowI)[j];

						newMat->addRow(newRow);
					}
				}
			}

			if (_rowNames != NULL) {
				newMat->_rowNames = new vector<string>();
				for (ullint i = 0; i < _rowNames->size(); ++i)
					if (deleteRows == NULL || deleteRows->find(i) == deleteRows->end())
						newMat->_rowNames->push_back((*_rowNames)[i]);
			}

			if (_colNames != NULL) {
				newMat->_colNames = new vector<string>();
				for (ullint j = 0; j < _colNames->size(); ++j)
					if (removeColumns->empty() || removeColumns->find(j) == removeColumns->end())
						newMat->_colNames->push_back((*_colNames)[j]);
			}

			delete removeColumns;

			return newMat;
		}

		NamedMatrix& normalize(bool normalizeMax = false) {
			RealType normFactor;
			if (normalizeMax)
				normFactor = max();
			else
				normFactor = sum();

			for (typename PARENT::const_iterator it = _mat->begin(); it != _mat->end(); ++it)
				**it /= normFactor;

			return *this;
		}

		inline NamedMatrix& operator=(const RealType& val) {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				**it = val;

			return *this;
		}

		inline NamedMatrix& operator+=(const RealType& add) {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				**it += add;

			return *this;
		}

		inline NamedMatrix& operator-=(const RealType& subtract) {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				**it -= subtract;

			return *this;
		}

		inline NamedMatrix& operator*=(const RealType& multiply) {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				**it *= multiply;

			return *this;
		}

		inline NamedMatrix& operator/=(const RealType& divide) {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				**it /= divide;

			return *this;
		}

		inline NamedMatrix& scaleByRowSums() {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				(*it)->scaleBySum();

			return *this;
		}

		inline NamedMatrix& scaleByColumnSums() {
			ullintPair dims = getDims();
			const ullint nrow = dims.first;
			const ullint ncol = dims.second;

			for (ullint j = 0; j < ncol; ++j) {
				DistributionStatistics<RealType> ds;
				for (ullint i = 0; i < nrow; ++i)
					ds.observeVal((*this)(i,j));
				RealType sum = ds.getSum();

				for (ullint i = 0; i < nrow; ++i) {
					RealType& val = (*this)(i,j);
					val /= sum;
				}
			}

			return *this;
		}

		inline NamedMatrix& log10(const RealType& pseudoCountValForLog10Input) {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				(*it)->log10(pseudoCountValForLog10Input);

			return *this;
		}

		inline NamedMatrix& centerRowsByMean(bool scaleByStdDev = false) {
			for (typename PARENT::iterator it = _mat->begin(); it != _mat->end(); ++it)
				(*it)->centerByMean(scaleByStdDev);

			return *this;
		}

		inline NamedMatrix& centerColumnsByMean(bool scaleByStdDev = false) {
			ullintPair dims = getDims();
			const ullint nrow = dims.first;
			const ullint ncol = dims.second;

			for (ullint j = 0; j < ncol; ++j) {
				DistributionStatistics<RealType> ds(scaleByStdDev);
				for (ullint i = 0; i < nrow; ++i)
					ds.observeVal((*this)(i,j));

				RealType mean = ds.getMean();
				RealType stdDev = ds.getStdDev();

				for (ullint i = 0; i < nrow; ++i) {
					RealType& val = (*this)(i,j);
					val -= mean;
					if (scaleByStdDev && stdDev != 0.0)
						val /= stdDev;
				}
			}

			return *this;
		}

		RealType sum() const {
			RealType sum = 0;
			for (typename PARENT::const_iterator it = _mat->begin(); it != _mat->end(); ++it)
				sum += (*it)->sum();

			return sum;
		}

		RealType max() const {
			typename PARENT::const_iterator it = _mat->begin();
			RealType max = (*it)->max();
			++it;

			for (; it != _mat->end(); ++it) {
				RealType curMax = (*it)->max();
				if (max < curMax)
					max = curMax;
			}

			return max;
		}

		RealType min() const {
			typename PARENT::const_iterator it = _mat->begin();
			RealType min = (*it)->min();
			++it;

			for (; it != _mat->end(); ++it) {
				RealType curMin = (*it)->min();
				if (curMin < min)
					min = curMin;
			}

			return min;
		}

		RealType columnSum(const ullint col) {
			RealType s = 0.0;
			for (ullint i = 0; i < nrow(); ++i)
				s += (*this)(i,col);

			return s;
		}

		void clear() {
			this->clearData();
		}

		typedef RealType (*MarginalizeFunc) (const RealType& v1, const RealType& v2);

		static RealType sum(const RealType& v1, const RealType& v2) {
			return v1 + v2;
		}

		static RealType max(const RealType& v1, const RealType& v2) {
			return (v1 > v2) ? v1 : v2;
		}

		static NamedMatrix<RealType>* multiplyMarginalizeMatrices(const NamedMatrix<RealType>& mat1, const NamedMatrix<RealType>& mat2, MarginalizeFunc margFunc = sum) {
			if (mat1.ncol() != mat2.nrow())
				throw new Exception("Cannot combine-marginalize matrices m1, m2 where ncol(m1) != nrow(m2)");
			const ullint commonDim = mat1.ncol();

			const ullint nrow = mat1.nrow();
			const ullint ncol = mat2.ncol();
			NamedMatrix<RealType>* combinedMat = new NamedMatrix<RealType>(nrow, ncol, 0.0);

			for (ullint i = 0; i < nrow; ++i) {
				for (ullint j = 0; j < ncol; ++j) {
					RealType val = 0;
					for (ullint k = 0; k < commonDim; ++k) {
						val = margFunc(val, mat1(i,k) * mat2(k,j));
					}
					(*combinedMat)(i,j) = val;
				}
			}

			return combinedMat;
		}

		NamedMatrix<RealType>* operator*(const NamedMatrix<RealType>& mat) const {
			return multiplyMarginalizeMatrices(*this, mat);
		}

		NamedMatrix<RealType>* operator+(const NamedMatrix<RealType>& mat) const {
			ullintPair dims = this->getDims();
			if (dims != mat.getDims())
				throw new Exception("Cannot add matrices of unequal dimensions");

			NamedMatrix<RealType>* sum = new NamedMatrix<RealType>(dims.first, dims.second);
			for (ullint i = 0; i < dims.first; ++i)
				for (ullint j = 0; j < dims.second; ++j)
					(*sum)(i,j) = (*this)(i,j) + mat(i,j);

			return sum;
		}

		NamedMatrix<RealType>* operator-(const NamedMatrix<RealType>& mat) const {
			ullintPair dims = this->getDims();
			if (dims != mat.getDims())
				throw new Exception("Cannot subtract matrices of unequal dimensions");

			NamedMatrix<RealType>* diff = new NamedMatrix<RealType>(dims.first, dims.second);
			for (ullint i = 0; i < dims.first; ++i)
				for (ullint j = 0; j < dims.second; ++j)
					(*diff)(i,j) = (*this)(i,j) - mat(i,j);

			return diff;
		}

		inline friend ostream& operator<<(ostream& stream, const NamedMatrix& nm) {
			return nm.printFixedWidth(stream, nm.namesAreOn());
		}

		virtual ostream& printDelimited(ostream& stream, bool printRowNames = true, bool printColNames = true, char delim = '\t') const {
			HMM_PP::ullintPair dims = getDims();
			const ullint nrow = dims.first;
			const ullint ncol = dims.second;

			if (printColNames) {
				stream << getMatrixName();
				for (ullint j = 0; j < ncol; ++j)
					stream << delim << colName(j);
				stream << endl;
			}

			for (ullint i = 0; i < nrow; ++i) {
				if (printRowNames)
					stream << rowName(i);
				stream << delim;
				(*_mat)[i]->printDelimited(stream, false, delim) << endl;
			}

			return stream;
		}

		virtual ostream& printFixedWidth(ostream& stream, bool showLabels = true, ullint fixedWidth = 5) const {
			ullintPair dims = getDims();
			const ullint nrow = dims.first;
			const ullint ncol = dims.second;

			if (nrow == 0) return stream;

			calcDefaultFixedWidth(fixedWidth, showLabels, stream.precision());

			if (showLabels) {
				if (fixedWidth > 0)
					stream << setw(fixedWidth);
				stream << '\t';
				for (ullint j = 0; j < ncol; ++j) {
					if (fixedWidth > 0)
						stream << setw(fixedWidth);
					stream << colName(j);
					if (j < ncol - 1)
						stream << '\t';
				}
				stream << endl;
			}

			for (ullint i = 0; i < nrow; ++i)  {
				if (showLabels) {
					if (fixedWidth > 0)
						stream << setw(fixedWidth);
					stream << rowName(i);
				}
				stream << '\t';
				(*_mat)[i]->printFixedWidth(stream, false, fixedWidth) << endl;
			}

			return stream;
		}

	protected:
		void calcDefaultFixedWidth(ullint& fixedWidth, bool showLabels, ullint precision) const{
			if (fixedWidth == 0)
				return;

			if (showLabels) {
				const ullint col = ncol();
				for (ullint j = 0; j < col; ++j) {
					string colJ = colName(j);
					if (colJ.size() > fixedWidth)
						fixedWidth = colJ.size();
				}
			}

			// TODO: a hacky estimate [5 characters: decimal point and 'e-10']:
			ullint requiredWith = precision + 5;
			if (fixedWidth < requiredWith)
				fixedWidth = requiredWith;
		}

	private:
		mutable PARENT* _mat; // _mat is mutable to permit lazy initialization of row vectors
		ullint _ncol;
		RealType _defaultValue;

		string _name;

		vector<string>* _rowNames;
		vector<string>* _colNames;
	};


	template<class RealType>
	const string NamedMatrix<RealType>::TRANSPOSE_SUFFIX = "_T";

	template<class RealType>
	NamedMatrix<RealType>* NamedVector<RealType>::diag() const {
		ullint n = this->size();
		NamedMatrix<RealType>* diag = new NamedMatrix<RealType>(n, n, 0.0);

		for (ullint i = 0; i < n; ++i)
			(*diag)(i,i) = (*this)[i];

		return diag;
	}

	template<class RealType>
	NamedMatrix<RealType>* NamedVector<RealType>::asMatrix(bool createColumnVector) const {
		NamedMatrix<RealType>* mat = new NamedMatrix<RealType>(*this);

		if (createColumnVector) {
			NamedMatrix<RealType>* tmpMat = mat->transpose();
			delete mat;
			mat = tmpMat;
		}

		return mat;
	}
}

#endif
