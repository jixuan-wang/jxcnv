#ifndef __DATA_OUPUTTER_H__
#define __DATA_OUPUTTER_H__

namespace HMM_PP {

	template<class DataType>
	class DataOutputter {

	public:
		DataOutputter() {}
		virtual ~DataOutputter() {}

		virtual void printDataType(const DataType* d) = 0;
	};
}

#endif
