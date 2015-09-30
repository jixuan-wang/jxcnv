#include "utils.hpp"
#include "Simplex.hpp"
#include "CRandom.hpp"
#include "Exception.hpp"
#include "params.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

const double HMM_PP::Simplex::EPS = 1e-8;

HMM_PP::Simplex::Simplex() {
	CRandom::srand(time(0));

	set_tolerance(0.00001);
	set_maxiter(2000);

	set_param(0);
	_func = NULL;
	_data = NULL;
	_iter = 0;
	_verbose = false;
}

void HMM_PP::Simplex::set_param(const uint n) {
	_ndim = n;
	_x.resize(n);
	_y.resize(n+1);
	_p.setDims(n+1, n);
	_psum.resize(n);
	_ptry.resize(n);
}

void HMM_PP::Simplex::set_start(const BaseRealVec& s) {
	// implicitly sets the # of parameters from 's' and
	// appropriately sizes all of the structures here
	set_param(s.size());

	for (uint i = 0; i < _ndim+1; i++)
		for (uint j = 0; j < _ndim; j++)
			_p(i,j) = s[j];

	randomize(0.5);
}

void HMM_PP::Simplex::randomize(const double fac) {
	for (uint i=0; i < _ndim+1; i++)
		for (uint j=0; j < _ndim; j++)
			_p(i,j) += (1 + fabs(_p(i,j))) * (CRandom::rand() - 0.5)* fac;
}

HMM_PP::BaseRealVec HMM_PP::Simplex::estimates() {
	// Take first row of simplex structure
	BaseRealVec parameters(_ndim);
	for (uint i = 0 ; i < _ndim ; i++)
		parameters[i] = _p(1,i);
	return parameters;
}

bool HMM_PP::Simplex::optimize(const uint rounds) {
	for (uint i = 0; i < rounds; i++) {
		if (i > 0) randomize(0.2);
		if (!optimize()) return false;
		cerr << "Round " << i << '\t' << _iter << ";\t";

		BaseRealVec e = estimates();
		for (uint j = 0; j < e.size(); j++)
			cerr << " " << e[j];
		cerr << endl;
		set_start(estimates());
	}

	return true;
}

bool HMM_PP::Simplex::optimize() {
	if (_func == NULL || _data == NULL) return false;

	// Initialize 'y':
	for (uint i = 0; i < _ndim+1; i++) {
		for (uint j = 0; j < _ndim; j++)
			_x[j] = _p(i,j);
		_y[i] = _func(_x, _data);
	}

	// Optimize:
	const double TINY = 1.0e-10;

	// p : n+1 rows, n columns
	// y : vector : n+1 elements
	uint ihi, ilo, inhi, mpts = _ndim+1;
	double rtol, ysave, ytry;
	_iter = 0;

	for (uint j=0; j < _ndim; j++)
		_psum[j] = _p.columnSum(j);

	for(;;) {
		ilo=0;
		ihi=_y[0]>_y[1] ? (inhi=1,0) : (inhi=0,1);

		for(uint i=0; i < mpts; i++) {
			if (_y[i] <= _y[ilo])
				ilo=i;

			if (_y[i] > _y[ihi]) {
				inhi=ihi;
				ihi=i;
			}
			else
				if (_y[i] > _y[inhi]&& i != ihi)
					inhi=i;
		}

		rtol = 2.0 * fabs(_y[ihi]-_y[ilo]) / (fabs(_y[ihi]) + fabs(_y[ilo]) + TINY);

		if (rtol < _ftol) {
			swap(_y[0], _y[ilo]);

			for (uint i=0; i < _ndim; i++)
				swap(_p(0,i), _p(ilo,i));
			break;
		}

		// Exceeded maximum number of iterations?
		if (_iter >= _NMAX) return false;

		_iter += 2;

		if (true || _verbose)  {
			cerr << "\n-2LL = " << _y[0] / 1000 << " (" << _iter << " iter) : param = \n";

			for (uint i=0; i < mpts; i++) {
				for (uint j=0; j < _ndim; j++)
					cerr << '\t' << _p(i,j);
				cerr << endl;
			}
			cerr << endl;
		}

		ytry = amotry(ihi, -1.0);

		if (ytry <= _y[ilo])
			ytry = amotry(ihi, 2.0);
		else if (ytry >= _y[inhi])  {
			ysave = _y[ihi];

			ytry = amotry(ihi, 0.5);

			if (ytry >= ysave)  {
				for (uint i=0; i < mpts; i++)  {
					if (i != ilo)  {
						//
						// TODO: Is the **DOUBLE** assignment below [ _p(i,j) = _psum[j] = ] correct?????
						//
						for(uint j=0; j < _ndim; j++)
							_p(i,j) = _psum[j] = 0.5 * (_p(i,j)+_p(ilo,j));
						_y[i]= _func(_psum, _data);
					}
				}
				_iter += _ndim;

				for (uint j=0; j < _ndim; j++)
					_psum[j] = _p.columnSum(j);
			}
		}
		else
			--_iter;
	}

	return true;
}

double HMM_PP::Simplex::amotry(const uint ihi, const double fac) {
	double fac1 = (1.0 - fac) / (double)_ndim;
	double fac2 = fac1 - fac;

	for(uint j = 0; j < _ndim; j++)
		_ptry[j] = _psum[j]* fac1 - _p(ihi,j)* fac2;

	double ytry = _func(_ptry, _data);

	if (!HMM_PP::realNum(ytry)) {
		cerr << "param:" << endl
				<< _ptry << endl;
		throw new Exception("function returned NaN in Nelder-Mead");
	}

	if (ytry < _y[ihi]) {
		_y[ihi]=ytry;

		for(uint j=0; j < _ndim; j++)  {
			_psum[j] += _ptry[j] - _p(ihi,j);
			_p(ihi,j) = _ptry[j];
		}
	}
	return ytry;
}
