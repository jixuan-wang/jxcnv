#ifndef __SIMPLEX_H__
#define __SIMPLEX_H__

#include "params.hpp"
#include "NamedVector.hpp"
#include "NamedMatrix.hpp"
#include "PreciseNonNegativeReal.hpp"

#include <cmath>
using namespace std;

namespace HMM_PP {

	class Simplex {
		typedef void data_t;
		typedef BaseReal (*func_t)(BaseRealVec&, data_t*);

	public:
		Simplex();

		// Initialize
		void set_function(func_t f) { _func = f; }

		void set_data(data_t* d) { _data = d ; }

		void set_start(const BaseRealVec&);

		// Drive
		bool optimize(const uint);

		bool optimize();

		BaseRealVec estimates();

		real likelihood() { return _y[0]; }

		// Misc.
		void set_tolerance(const double tol) { _ftol = tol; }

		void set_maxiter(const uint n) { _NMAX = n; }

		uint n_param() const { return _ndim; }

		uint n_iter() const { return _iter; }

		void randomize(const double fac = 5);

	private:
		static const double EPS;

		double _ftol;
		uint _ndim;
		uint _iter;
		uint _NMAX;
		bool _verbose;

		func_t _func;
		data_t* _data;

		BaseRealMat _p;
		BaseRealVec _x;
		BaseRealVec _y;
		BaseRealVec _psum;
		BaseRealVec _ptry;

		void set_param(const uint n);

		double amotry(const uint ihi, const double fac);
	};
}

#endif

