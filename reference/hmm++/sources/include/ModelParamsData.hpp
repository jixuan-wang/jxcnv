#ifndef __MODEL_PARAMS_DATA_H__
#define __MODEL_PARAMS_DATA_H__

#include "ModelParams.hpp"
#include "DataLoader.hpp"

namespace HMM_PP {

	template<class Type>
	class ModelParamsData : virtual public ModelParams {

	public:
		typedef Type DataType;

		ModelParamsData(HiddenStateParams* params, DataLoader<Type>* dataLoader) : ModelParams(params), _dataLoader(dataLoader) {}
		virtual ~ModelParamsData() = 0;

		virtual Type* getNewData() const {
			return new Type();
		}

		virtual bool hasNextLoadData() const {
			return _dataLoader != NULL && _dataLoader->hasNext();
		}

		virtual Type* loadNextData() const {
			if (!hasNextLoadData())
				throw new Exception("Cannot load from DataLoader when !hasNextLoadData()");

			return _dataLoader->next();
		}

	protected:
		static Type* castData(Data* d) {
			return dynamic_cast<Type*>(d);
		}

		static const Type* castData(const Data* d) {
			return dynamic_cast<const Type*>(d);
		}

		DataLoader<Type>* _dataLoader;
	};

	template<class Type>
	ModelParamsData<Type>::~ModelParamsData() {}
}

#endif
