#ifndef __UNIVARIATE_DATA_H__
#define __UNIVARIATE_DATA_H__

#include "params.hpp"
#include "Data.hpp"
#include "NamedVector.hpp"

namespace HMM_PP {

	// Observations for a particular sequence
	template<class BaseType>
	class UnivariateData : virtual public Data {

	public:
		UnivariateData();
		virtual ~UnivariateData();

	protected:
		virtual void clearObservations();

	public:
		virtual void setNumObservations(const uint nobs);

		void setDatapoint(const uint i, const BaseType& t);
		void addDatapoint(const BaseType& t);

		// get datapoint at time 't'
		BaseType val(const uint t) const;

		virtual void printDatapoint(ostream& stream, const uint t) const;

		const NamedVector<BaseType>& getValues() const { return *_data; }

		NamedVector<BaseType>& getNonConstValues() { return *_data; }
		NamedVector<BaseType>* getValuesTransferOwner() { NamedVector<BaseType>* tmp = _data; _data = NULL; return tmp; }

	private:
		NamedVector<BaseType>* _data;
	};

	template<class BaseType>
	UnivariateData<BaseType>::UnivariateData()
	: Data(), _data(new NamedVector<BaseType>()) {
	}

	template<class BaseType>
	UnivariateData<BaseType>::~UnivariateData() {
		if (_data != NULL)
			delete _data;
	}

	template<class BaseType>
	void UnivariateData<BaseType>::clearObservations() {
		_data->clear();
	}

	template<class BaseType>
	void UnivariateData<BaseType>::setNumObservations(const uint nobs) {
		Data::setNumObservations(nobs);
		_data->clear();
		_data->resize(_n);
	}

	template<class BaseType>
	void UnivariateData<BaseType>::setDatapoint(const uint i, const BaseType& t) {
		if (i < _n)
			(*_data)[i] = t;
	}

	template<class BaseType>
	void UnivariateData<BaseType>::addDatapoint(const BaseType& t) {
		_data->addElement(t);
		++_n;
	}

	template<class BaseType>
	BaseType UnivariateData<BaseType>::val(const uint t) const {
		return (*_data)[t];
	}

	template<class BaseType>
	void UnivariateData<BaseType>::printDatapoint(ostream& stream, const uint t) const {
		stream << val(t);
	}
}

#endif
