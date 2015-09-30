#ifndef __PARAMS_H__
#define __PARAMS_H__

#include <cmath>
#include <climits>
using namespace std;


typedef unsigned int uint;
#define UINT_INFINITY UINT_MAX

typedef unsigned long long int ullint;
#define ULLINT_INFINITY ULLONG_MAX

typedef long int lint;
#define LINT_INFINITY LONG_MAX


typedef long double BaseReal;


namespace HMM_PP {
	template<class RealType> class NamedVector;
	template<class RealType> class NamedMatrix;

	typedef NamedVector<double> DoubleVec;
	typedef NamedMatrix<double> DoubleMat;

	typedef NamedVector<BaseReal> BaseRealVec;
	typedef NamedMatrix<BaseReal> BaseRealMat;
}


#endif
