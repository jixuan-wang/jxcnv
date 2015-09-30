#ifndef __NAMED_VECTOR_H__
#define __NAMED_VECTOR_H__

#include "params.hpp"
#include "Exception.hpp"
#include "DistributionStatistics.hpp"

#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <utility>
#include <cmath>
using namespace std;

namespace HMM_PP {
	template<class RealType> class NamedMatrix;

	template<class RealType>
	class NamedVector {

	public:
		typedef vector<RealType> PARENT;

		explicit NamedVector() : _vec(new PARENT()), _constLength(false), _names(NULL) {}
		explicit NamedVector(const ullint s, const RealType val = RealType(), bool constLength = false) : _vec(new PARENT(s, val)), _constLength(constLength), _names(NULL) {}
		NamedVector(const PARENT& copy, bool constLength = false) : _vec(new PARENT(copy)), _constLength(constLength), _names(NULL) {}
		NamedVector(const NamedVector& copy, bool constLength = false) : _vec(new PARENT(*(copy._vec))), _constLength(constLength), _names(NULL) {}

		virtual ~NamedVector() {
			if (_vec != NULL)
				delete _vec;

			if (_names != NULL)
				delete _names;
		}

		void setConstantLength() {
			_constLength = true;
		}

		NamedVector& operator=(const PARENT& v) {
			if (_constLength && v.size() != this->size())
				throw new Exception("Cannot modify the length of constant-length NamedVector");

			*_vec = v;
			return *this;
		}

		NamedVector& operator=(const NamedVector& v) {
			return operator=(*(v._vec));
		}

		void addElement(const RealType& v) {
			if (_constLength)
				throw new Exception("Cannot modify the length of constant-length NamedVector");

			_vec->push_back(v);
		}

		void resize(const ullint s, const RealType val = RealType()) {
			if (_constLength && s != _vec->size())
				throw new Exception("Cannot modify the length of constant-length NamedVector");

			_vec->resize(s, val);
		}

		inline ullint size() const { return _vec->size(); }
		inline bool empty() const { return _vec->empty(); }

		const RealType& operator[](const ullint i) const { return (*_vec)[i]; }
		RealType& operator[](const ullint i) { return (*_vec)[i]; }

	private:
		void turnOnNames() {
			ullint sz = _vec->size();

			if (sz > 0) {
				if (_names == NULL)
					_names = new vector<string>(sz);
				else
					_names->resize(sz);
			}
			else if (_names == NULL)
				_names = new vector<string>();
		}

		bool namesAreOn() const {
			return _names != NULL;
		}

	public:
#define CHECK_INDEX \
		if (i > _vec->size()) { \
			stringstream str; \
			str << "Index " << i << " is out of bounds"; \
			throw new Exception(str.str()); \
		}

		void setName(const ullint i, const string& s) {
			CHECK_INDEX

			if (!namesAreOn())
				turnOnNames();

			ullint minSz = i + 1;
			if (_names->size() < minSz)
				_names->resize(_vec->size());

			(*_names)[i] = s;
		}

		string name(const ullint i) const {
			CHECK_INDEX

			ullint minSz = i + 1;
			if (_names == NULL || _names->size() < minSz) {
				stringstream str;
				str << "[" << minSz << "]";
				return str.str();
			}

			return (*_names)[i];
		}

		NamedVector& normalize(bool normalizeMax = false) {
			RealType normFactor;
			if (normalizeMax)
				normFactor = max();
			else
				normFactor = sum();

			*this /= normFactor;

			return *this;
		}

		inline NamedVector& operator=(const RealType& val) {
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it)
				*it = val;

			return *this;
		}

		inline NamedVector& operator+=(const RealType& add) {
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it)
				*it += add;

			return *this;
		}

		inline NamedVector& operator-=(const RealType& subtract) {
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it)
				*it -= subtract;

			return *this;
		}

		inline NamedVector& operator*=(const RealType& multiply) {
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it)
				*it *= multiply;

			return *this;
		}

		inline NamedVector& operator/=(const RealType& divide) {
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it)
				*it /= divide;

			return *this;
		}

		NamedVector& operator+=(const PARENT& add) {
			if (add.size() != this->size())
				throw new Exception("Cannot calculate the dot-addition of vectors of differing lengths!");

			typename PARENT::const_iterator addIt = add.begin();
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it) {
				*it += *addIt;
				++addIt;
			}

			return *this;
		}

		NamedVector& operator+=(const NamedVector& v) {
			return operator+=(*(v._vec));
		}

		NamedVector& operator-=(const PARENT& subtract) {
			if (subtract.size() != this->size())
				throw new Exception("Cannot calculate the dot-subtraction of vectors of differing lengths!");

			typename PARENT::const_iterator subtractIt = subtract.begin();
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it) {
				*it -= *subtractIt;
				++subtractIt;
			}

			return *this;
		}

		NamedVector& operator-=(const NamedVector& v) {
			return operator-=(*(v._vec));
		}

		NamedVector& operator*=(const PARENT& multiply) {
			if (multiply.size() != this->size())
				throw new Exception("Cannot calculate the dot-multiplication of vectors of differing lengths!");

			typename PARENT::const_iterator multIt = multiply.begin();
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it) {
				*it *= *multIt;
				++multIt;
			}

			return *this;
		}

		NamedVector& operator*=(const NamedVector& v) {
			return operator*=(*(v._vec));
		}

		NamedVector& operator/=(const PARENT& divide) {
			if (divide.size() != this->size())
				throw new Exception("Cannot calculate the dot-division of vectors of differing lengths!");

			typename PARENT::const_iterator divideIt = divide.begin();
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it) {
				*it /= *divideIt;
				++divideIt;
			}

			return *this;
		}

		NamedVector& operator/=(const NamedVector& v) {
			return operator/=(*(v._vec));
		}

		inline NamedVector& scaleBySum() {
			DistributionStatistics<RealType> ds;

			for (typename PARENT::const_iterator it = _vec->begin(); it != _vec->end(); ++it)
				ds.observeVal(*it);
			RealType sum = ds.getSum();

			return (*this /= sum);
		}

		NamedVector& log10(const RealType& pseudoCountValForLog10Input) {
			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it) {
				RealType val = *it;
				if (val <= 0)
					val = 0;
				val += pseudoCountValForLog10Input;

				*it = std::log10(val);
			}

			return *this;
		}

		inline NamedVector& centerByMean(bool scaleByStdDev = false) {
			DistributionStatistics<RealType> ds(scaleByStdDev);

			for (typename PARENT::const_iterator it = _vec->begin(); it != _vec->end(); ++it)
				ds.observeVal(*it);

			RealType mean = ds.getMean();
			RealType stdDev = ds.getStdDev();

			for (typename PARENT::iterator it = _vec->begin(); it != _vec->end(); ++it) {
				*it -= mean;
				if (scaleByStdDev && stdDev != 0.0)
					*it /= stdDev;
			}

			return *this;
		}

		RealType sum() const {
			RealType sum = 0;
			for (typename PARENT::const_iterator it = _vec->begin(); it != _vec->end(); ++it)
				sum += *it;

			return sum;
		}

		pair<RealType, ullint> calcMaxArgmax() const {
			ullint sz = this->size();
			if (sz <= 0)
				throw new Exception("Cannot calculate maximum of empty vector!");

			ullint argmax = 0;
			RealType max = (*this)[argmax];
			for (ullint i = 0; i < sz; ++i) {
				RealType val = (*this)[i];
				if (val > max) {
					max = val;
					argmax = i;
				}
			}

			return pair<RealType, ullint>(max, argmax);
		}

		RealType max() const {
			return calcMaxArgmax().first;
		}

		pair<RealType, ullint> calcMinArgmin() const {
			ullint sz = this->size();
			if (sz <= 0)
				throw new Exception("Cannot calculate minimum of empty vector!");

			ullint argmin = 0;
			RealType min = (*this)[argmin];
			for (ullint i = 0; i < sz; ++i) {
				RealType val = (*this)[i];
				if (val < min) {
					min = val;
					argmin = i;
				}
			}

			return pair<RealType, ullint>(min, argmin);
		}

		RealType min() const {
			return calcMinArgmin().first;
		}

		NamedMatrix<RealType>* diag() const;
		NamedMatrix<RealType>* asMatrix(bool createColumnVector = false) const;

		void clear() {
			_vec->clear();
		}

		inline friend ostream& operator<<(ostream& stream, const NamedVector& nv) {
			return nv.printFixedWidth(stream, nv.namesAreOn());
		}

		virtual ostream& printDelimited(ostream& stream, bool printNames = true, char delim = '\t') const {
			ullint sz = this->size();

			if (printNames) {
				for (ullint i = 0; i < sz; ++i) {
					stream << name(i);
					if (i < sz - 1)
						stream << delim;
				}
				stream << endl;
			}

			for (ullint i = 0; i < sz; ++i) {
				stream << (*this)[i];
				if (i < sz - 1)
					stream << delim;
			}

			return stream;
		}

		virtual ostream& printFixedWidth(ostream& stream, bool showLabels = true, const int fixedWidth = 0) const {
			ullint sz = this->size();

			if (showLabels) {
				for (ullint i = 0; i < sz; i++) {
					if (fixedWidth > 0)
						stream << setw(fixedWidth);
					stream << name(i);
					if (i < sz - 1)
						stream << '\t';
				}
				stream << endl;
			}

			for (ullint i = 0; i < sz; i++) {
				if (fixedWidth > 0)
					stream << setw(fixedWidth);
				stream << (*this)[i];
				if (i < sz - 1)
					stream << '\t';
			}

			return stream;
		}

	private:
		PARENT* _vec;
		bool _constLength;

		vector<string>* _names;
	};
}

#endif
