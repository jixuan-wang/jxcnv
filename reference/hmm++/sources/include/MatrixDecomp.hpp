#ifndef __MATRIX_DECOMP_H__
#define __MATRIX_DECOMP_H__

#include "NamedVector.hpp"
#include "NamedMatrix.hpp"
#include "LAPACKvector.hpp"

#include <set>
using namespace std;

namespace HMM_PP {
	class MatrixDecomp {
	public:
		template<class RealType>
		class SVDecomp {
		public:
			SVDecomp(DoubleVec* singularVals, NamedMatrix<RealType>* leftSingVectors, NamedMatrix<RealType>* rightSingVectors);
			~SVDecomp();
			void deleteData();

			DoubleVec* D; // singularVals
			NamedMatrix<RealType>* U; // leftSingVectors
			NamedMatrix<RealType>* V; // rightSingVectors
		};

		template<class RealType>
		class EigenDecomp {
		public:
			EigenDecomp(DoubleVec* eigenValsReal, DoubleVec* eigenValsImaginary, NamedMatrix<RealType>* eigenVectors);
			~EigenDecomp();
			void deleteData();

			DoubleVec* eigenValsReal;
			DoubleVec* eigenValsImaginary;
			NamedMatrix<RealType>* eigenVectors;
		};

		/*
		 * Computes the SVD decomposition of matrix X:
		 * X = U * D * V'
		 *
		 * where U and V are orthogonal, V' means _V transposed_, and D is a
		 * diagonal matrix with the singular values D[i,i].  Equivalently, D = U' X V.
		 *
		 * The SVDecomp object returned contains:
		 * D: a vector containing the singular values of X, sorted in decreasing order
		 * U: a matrix whose columns contain the corresponding left singular vectors of X
		 * V: a matrix whose columns contain the corresponding right singular vectors of X
		 */
		template<class RealType>
		static SVDecomp<RealType> svd(NamedMatrix<RealType>* X, bool returnVtranspose = false, string workDir = "", bool deleteX = true);

		/*
		 * Computes the Right-Eigen decomposition of square matrix X.
		 *
		 * The EigenDecomp object returned contains:
		 * eigenValsReal: a vector containing the *real* parts of the eigenvalues of X
		 * eigenValsImaginary: a vector containing the *imaginary* parts of the corresponding eigenvalues of X
		 * eigenVectors: a matrix whose columns contain the corresponding right eigenvectors of X
		 */
		template<class RealType>
		static EigenDecomp<RealType> rightEigen(NamedMatrix<RealType>* X, string workDir = "", bool deleteX = true);

		static void initSignalHandlers();

	private:
		template<class RealType>
		static LAPACKvector<double>* matrixToLAPACKvector(NamedMatrix<RealType>* mat, bool deleteMat, const string& workDir, bool transpose = false);

		template<class RealType>
		static NamedMatrix<RealType>* LAPACKvectorToMatrix(const double* matVec, ullint m, ullint n, bool transposeMatrix = false);

		static LAPACKvector<double>* allocateVector(ullint size, const string workDir = "");
		static void deleteAllocatedVector(HMM_PP::LAPACKvector<double>* vec);

		static void freeMappedMemory(int sig);
		static set<LAPACKvector<double>*> _mappedMem;
	};




	/*
	 * NOTE:
	 * As described here:
	 * http://bytes.com/topic/c/answers/453169-dynamic-arrays-convert-vector-array
	 *
	 * NOTE: We can pass &vec[0] as a pointer to a contiguous memory block that contains the values in the array underlying std::vector vec
	 */
	extern "C" int dgesdd_(const char* jobz, lint* m, lint* n, double* a,
			lint* lda, double* s, double* u, lint* ldu,
			double* vt, lint* ldvt, double* work, lint* lwork,
			lint* iwork, lint* info);

	// "S": the first min(M,N) columns of U and the first min(M,N) rows of V**T are returned in the arrays U and VT:
	#define LAPACKsvd dgesdd_("S", &m, &n, Xvec, &m, &(*D)[0], Uvec, &m, VtransposeVec->rawData(), &m, work, &lwork, iwork, &info)

	template<class RealType>
	HMM_PP::MatrixDecomp::SVDecomp<RealType> HMM_PP::MatrixDecomp::svd(NamedMatrix<RealType>* X, bool returnVtranspose, string workDir, bool deleteX) {
		ullintPair m_n = X->getDims();
		lint m = m_n.first;
		lint n = m_n.second;
		if (m <= 0 || n <= 0)
			throw new Exception("Cannot SVDecompose an empty matrix");

		bool swapped = false;
		LAPACKvector<double>* lVec = NULL;

		if (m <= n) {
			lVec = matrixToLAPACKvector(X, deleteX, workDir);
		}
		else {
			swap(m, n);
			swapped = true;

			lVec = matrixToLAPACKvector(X, deleteX, workDir, true);
		}
		double* Xvec = lVec->rawData();

		DoubleVec* D = new DoubleVec(m);
		double* Uvec = new double[m * m];
		LAPACKvector<double>* VtransposeVec = allocateVector(m * n, workDir);

		double* work = new double[1];
		lint lwork = -1;
		lint* iwork = new lint[8 * m];
		lint info = 0;

		// Determine size of workspace needed for SVD:
		LAPACKsvd;
		if (info != 0)
			throw new Exception("Error running FIRST dgesdd_");

		lwork = static_cast<lint>(work[0]);
		delete[] work;
		LAPACKvector<double>* workVec = allocateVector(lwork, workDir);
		work = workVec->rawData();

		// Perform actual SVD:
		LAPACKsvd;
		if (info != 0)
			throw new Exception("Error running SECOND dgesdd_");

		delete[] iwork;
		deleteAllocatedVector(workVec);
		deleteAllocatedVector(lVec);

		NamedMatrix<RealType>* U = LAPACKvectorToMatrix<RealType>(Uvec, m, m);
		U->setMatrixName("U");
		delete[] Uvec;

		NamedMatrix<RealType>* V = LAPACKvectorToMatrix<RealType>(VtransposeVec->rawData(), m, n, !returnVtranspose);
		string matName = "V";
		if (returnVtranspose)
			matName += NamedMatrix<RealType>::TRANSPOSE_SUFFIX;
		V->setMatrixName(matName);
		deleteAllocatedVector(VtransposeVec);

		if (swapped)
			swap(U, V);

		return SVDecomp<RealType>(D, U, V);
	}

	template<class RealType>
	HMM_PP::MatrixDecomp::SVDecomp<RealType>::SVDecomp(DoubleVec* singularVals, NamedMatrix<RealType>* leftSingVectors, NamedMatrix<RealType>* rightSingVectors)
	: D(singularVals), U(leftSingVectors), V(rightSingVectors) {
	}

	template<class RealType>
	HMM_PP::MatrixDecomp::SVDecomp<RealType>::~SVDecomp() {
	}

	template<class RealType>
	void HMM_PP::MatrixDecomp::SVDecomp<RealType>::deleteData() {
		if (D != NULL) {
			delete D;
			D = NULL;
		}

		if (U != NULL) {
			delete U;
			U = NULL;
		}

		if (V != NULL) {
			delete V;
			V = NULL;
		}
	}

	extern "C" int dgeev_(const char*, const char*, lint*,
			double*, lint*, double*, double*,
			double*, lint*, double*, lint*,
			double*, lint*, lint*);

	#define LAPACKeigen \
			dgeev_("N",      /* left eigenvectors are not computed */ \
			"V",      /* right eigenvectors are computed */ \
			&n ,      /* order of matrix */ \
			Xvec,    /* input matrix */ \
			&n,       /* LDA */ \
			&(*eigenValsReal)[0], /* real parts of the computed eigenvalues */ \
			&(*eigenValsImaginary)[0], /* imaginary parts of the computed eigenvalues */ \
			&DUMMY_VL, /* left eigenvectors */ \
			&DUMMY_LDVL, /* LDVL */ \
			eigenVecsVector, /* right eigenvectors */ \
			&n, /* LDVR */ \
			work,  /* Workspace */ \
			&lwork, /* size of workspace */ \
			&info)

	template<class RealType>
	HMM_PP::MatrixDecomp::EigenDecomp<RealType> HMM_PP::MatrixDecomp::rightEigen(NamedMatrix<RealType>* X, string workDir, bool deleteX) {
		ullintPair m_n = X->getDims();
		lint m = m_n.first;
		lint n = m_n.second;
		if (m != n || n <= 0)
			throw new Exception("Cannot Eigen-decompose a non-square or empty matrix");

		LAPACKvector<double>* lVec = matrixToLAPACKvector(X, deleteX, workDir);
		double* Xvec = lVec->rawData();

		DoubleVec* eigenValsReal = new DoubleVec(n, 0);
		DoubleVec* eigenValsImaginary = new DoubleVec(n, 0);

		double* eigenVecsVector = new double[n * n];

		double* work = new double[1];
		lint lwork = -1;
		lint info = 0;

		double DUMMY_VL;
		lint DUMMY_LDVL = 1;

		// Get workspace:
		LAPACKeigen;
		if (info != 0)
			throw new Exception("Error running FIRST dgeev_");

		// Assign workspace:
		lwork = static_cast<lint>(work[0]);
		delete[] work;
		work = new double[lwork];

		LAPACKeigen;
		if (info != 0)
			throw new Exception("Error running SECOND dgeev_");

		delete lVec;
		delete[] work;

		NamedMatrix<RealType>* eigenVecsMatrix = LAPACKvectorToMatrix<RealType>(eigenVecsVector, n, n);
		delete[] eigenVecsVector;

		return EigenDecomp<RealType>(eigenValsReal, eigenValsImaginary, eigenVecsMatrix);
	}

	template<class RealType>
	HMM_PP::MatrixDecomp::EigenDecomp<RealType>::EigenDecomp(DoubleVec* eigenValsReal, DoubleVec* eigenValsImaginary, NamedMatrix<RealType>* eigenVectors)
	: eigenValsReal(eigenValsReal), eigenValsImaginary(eigenValsImaginary), eigenVectors(eigenVectors) {
	}

	template<class RealType>
	HMM_PP::MatrixDecomp::EigenDecomp<RealType>::~EigenDecomp() {
	}

	template<class RealType>
	void HMM_PP::MatrixDecomp::EigenDecomp<RealType>::deleteData() {
		if (eigenValsReal != NULL) {
			delete eigenValsReal;
			eigenValsReal = NULL;
		}

		if (eigenValsImaginary != NULL) {
			delete eigenValsImaginary;
			eigenValsImaginary = NULL;
		}

		if (eigenVectors != NULL) {
			delete eigenVectors;
			eigenVectors = NULL;
		}
	}

	template<class RealType>
	HMM_PP::LAPACKvector<double>* HMM_PP::MatrixDecomp::matrixToLAPACKvector(NamedMatrix<RealType>* mat, bool deleteMat, const string& workDir, bool transpose) {
		HMM_PP::ullintPair m_n = mat->getDims();
		ullint m = m_n.first;  // # of rows
		ullint n = m_n.second; // # of columns

		HMM_PP::LAPACKvector<double>* retVec = allocateVector(m * n);
		double* matVec = retVec->rawData();

		ullint ind = 0;
		ullint row = m-1;
		while (true) {
			for (ullint col = 0; col < n; ++col) {
				if (!transpose) // LAPACK requires conversion to vector in column-major order:
					ind = col * m + row;
				else // want to transpose, so insert data in row-major order:
					ind = row * n + col;

				matVec[ind] = (*mat)(row, col);
			}

			// Delete as we go to save space:
			if (deleteMat)
				mat->dropLastRow();

			if (row == 0)
				break;
			else
				--row;
		}

		if (deleteMat)
			delete mat;

		return retVec;
	}

	template<class RealType>
	HMM_PP::NamedMatrix<RealType>* HMM_PP::MatrixDecomp::LAPACKvectorToMatrix(const double* matVec, ullint m, ullint n, bool transposeMatrix) {
		NamedMatrix<RealType>* mat = new NamedMatrix<RealType>();

		ullint k = 0;

		if (!transposeMatrix) { // use the default LAPACK column-major order
			mat->setDims(m, n);
			for (ullint col = 0; col < n; ++col)
				for (ullint row = 0; row < m; ++row)
					(*mat)(row, col) = matVec[k++];
		}
		else { // transpose the matrix, so insert in row-major order
			mat->setDims(n, m);
			for (ullint row = 0; row < n; ++row)
				for (ullint col = 0; col < m; ++col)
					(*mat)(row, col) = matVec[k++];
		}

		return mat;
	}


}

#endif
