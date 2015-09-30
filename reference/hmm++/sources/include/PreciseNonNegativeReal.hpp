#ifndef __PRECISE_NON_NEGATIVE_REAL_H__
#define __PRECISE_NON_NEGATIVE_REAL_H__

#include "Exception.hpp"
#include "params.hpp"

#include <cmath>
#include <sstream>
using namespace std;

#define real HMM_PP::PreciseNonNegativeReal<BaseReal>

#define RealVec NamedVector<real>
#define RealMat NamedMatrix<real>

/* PreciseNonNegativeReal permits arithmetic operations on NON-NEGATIVE real values
   with precision (prevents underflow by representing in log10 space).

   Note that cmath has std::log10() and std::pow() versions for various RealType data types
 */
namespace HMM_PP {

	template<class RealType>
	class PreciseNonNegativeReal {

	private:
		static const RealType REAL_ZERO;
		static const RealType REAL_ONE;
		static const RealType REAL_TEN;

	public:
		static const RealType REAL_EPSILON;

		static const RealType REAL_INFINITY;

		PreciseNonNegativeReal(RealType d = REAL_ZERO, bool isLog10 = false) {
			if (isLog10) {
				this->log10Value = d;
			}
			else {
				if (d < REAL_ZERO) {
					stringstream str;
					str << "non-log PreciseNonNegativeReal argument must be non-negative: " << d;
					throw new Exception(str.str());
				}
				this->log10Value = std::log10(d);
			}
		}

		PreciseNonNegativeReal(const PreciseNonNegativeReal& pd) {
			this->log10Value = pd.log10Value;
		}

		inline RealType getValue() const {
			return std::pow(REAL_TEN, log10Value);
		}

		inline RealType getLog10Value() const {
			return log10Value;
		}

		PreciseNonNegativeReal& operator=(const PreciseNonNegativeReal& other) {
			log10Value = other.log10Value;
			return *this;
		}

		PreciseNonNegativeReal operator+(const PreciseNonNegativeReal& other) const {
			PreciseNonNegativeReal tmp(*this);
			tmp += other;
			return tmp;
		}

		PreciseNonNegativeReal operator-(const PreciseNonNegativeReal& other) const {
			PreciseNonNegativeReal tmp(*this);
			tmp -= other;
			return tmp;
		}

		PreciseNonNegativeReal operator*(const PreciseNonNegativeReal& other) const {
			PreciseNonNegativeReal tmp(*this);
			tmp *= other;
			return tmp;
		}

		PreciseNonNegativeReal operator/(const PreciseNonNegativeReal& other) const {
			PreciseNonNegativeReal tmp(*this);
			tmp /= other;
			return tmp;
		}

		inline bool operator==(const PreciseNonNegativeReal& other) const {
			return (this->compareTo(other) == 0);
		}

		inline bool operator>(const PreciseNonNegativeReal& other) const {
			return (this->compareTo(other) > 0);
		}

		inline bool operator<(const PreciseNonNegativeReal& other) const {
			return (this->compareTo(other) < 0);
		}

		inline bool operator>=(const PreciseNonNegativeReal& other) const {
			return (this->compareTo(other) >= 0);
		}

		inline bool operator<=(const PreciseNonNegativeReal& other) const {
			return (this->compareTo(other) <= 0);
		}

		inline PreciseNonNegativeReal& operator+=(const PreciseNonNegativeReal& other) {
			log10Value = addInLogSpace(log10Value, other.log10Value);
			return *this;
		}

		inline PreciseNonNegativeReal& operator-=(const PreciseNonNegativeReal& other) {
			log10Value = absSubLog(log10Value, other.log10Value);
			return *this;
		}

		inline PreciseNonNegativeReal& operator*=(const PreciseNonNegativeReal& other) {
			log10Value += other.log10Value;
			return *this;
		}

		inline PreciseNonNegativeReal& operator/=(const PreciseNonNegativeReal& other) {
			log10Value -= other.log10Value;
			return *this;
		}

	private:
		RealType log10Value;

		// If x = log(a), y = log(b), returns log(a+b)
		static RealType addInLogSpace(RealType x, RealType y) {
			if (x == REAL_INFINITY || y == REAL_INFINITY)
				return REAL_INFINITY; // log(e^INFINITY + e^y) = INFINITY

			if (x == -REAL_INFINITY)
				return y;
			if (y == -REAL_INFINITY)
				return x;

			RealType maxVal, negDiff;
			if (x > y) {
				maxVal = x;
				negDiff = y - x;
			}
			else { // x <= y
				maxVal = y;
				negDiff = x - y;
			}

			// x + log(1+e^(y-x)) = log(a) + log(1+e^(log(b)-log(a))) = log(a) + log(1+b/a) = log(a+b)
			return maxVal + std::log10(REAL_ONE + std::pow(REAL_TEN, negDiff));
		}

		// If x = log(a), y = log(b), returns log |a-b|
		static RealType absSubLog(RealType x, RealType y) {
			if (x == -REAL_INFINITY && y == -REAL_INFINITY)
				return -REAL_INFINITY; // log |e^-INFINITY - e^-INFINITY| = log |0-0| = log(0) = -INFINITY

			RealType maxVal, negDiff;
			if (x > y) {
				maxVal = x;
				negDiff = y - x;
			}
			else { // x <= y
				maxVal = y;
				negDiff = x - y;
			}

			// x + log(1-e^(y-x)) = log(a) + log(1-e^(log(b)-log(a))) = log(a) + log(1-b/a) = a - b = |a-b|, since x >= y
			return maxVal + std::log10(REAL_ONE - std::pow(REAL_TEN, negDiff));
		}

		int compareTo(const PreciseNonNegativeReal& other) const {
			// Since log is monotonic: e^a R e^b <=> a R b, where R is one of: >, <, ==
			return compareReal(this->log10Value, other.log10Value);
		}

	public:
		static int compareReal(const RealType& r1, const RealType& r2) {
			RealType diff = r1 - r2;
			if (abs(diff) <= REAL_EPSILON)
				return 0; // r1 "==" r2

			if (diff > REAL_ZERO)
				return 1;
			else
				return -1;
		}

		static bool equalsReal(const RealType& r1, const RealType& r2) {
			return compareReal(r1, r2) == 0;
		}
	};

	template <class RealType>
	ostream& operator<<(ostream& stream, const PreciseNonNegativeReal<RealType>& pnnr) {
		stream << pnnr.getValue();
		return stream;
	}
}


template<class RealType>
const RealType HMM_PP::PreciseNonNegativeReal<RealType>::REAL_ZERO = 0.0;

template<class RealType>
const RealType HMM_PP::PreciseNonNegativeReal<RealType>::REAL_ONE = 1.0;

template<class RealType>
const RealType HMM_PP::PreciseNonNegativeReal<RealType>::REAL_TEN = 10.0;

template<class RealType>
const RealType HMM_PP::PreciseNonNegativeReal<RealType>::REAL_EPSILON = 1e-6;

template<class RealType>
const RealType HMM_PP::PreciseNonNegativeReal<RealType>::REAL_INFINITY = - std::log10(PreciseNonNegativeReal<RealType>::REAL_ZERO);


#endif
