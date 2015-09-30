#ifndef __LAPACK_VECTOR_H__
#define __LAPACK_VECTOR_H__

namespace HMM_PP {

	template<class dataType>
	class LAPACKvector {

	protected:
		LAPACKvector() {}

	public:
		virtual ~LAPACKvector() {}

		virtual const dataType& operator[](const ullint i) const = 0;
		virtual dataType& operator[](const ullint i) = 0;

		virtual dataType* rawData() = 0;
		virtual const dataType* rawData() const = 0;
	};


	template<class dataType>
	class StandardVector : public LAPACKvector<dataType> {
	public:
		StandardVector(ullint size) : _vec(new dataType[size]) {}
		virtual ~StandardVector() { delete[] _vec; }

		virtual const dataType& operator[](const ullint i) const { return _vec[i]; }
		virtual dataType& operator[](const ullint i) { return _vec[i]; }

		virtual dataType* rawData() { return _vec; }
		virtual const dataType* rawData() const { return _vec; }

	private:
		dataType* _vec;
	};
}

#endif
