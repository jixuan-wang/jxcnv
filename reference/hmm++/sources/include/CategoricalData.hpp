#ifndef __CATEGORICAL_DATA_H__
#define __CATEGORICAL_DATA_H__

#include "params.hpp"
#include "Data.hpp"

#include <map>
#include <string>
using namespace std;

namespace HMM_PP {
	class Model;

	// Observations for a particular sequence
	class CategoricalData : virtual public Data {

	public:
		CategoricalData() : Data() {}
		virtual ~CategoricalData() = 0;

		static uint ensureCategoryToType(const string& c);
		static uint categoryToType(const string& c);
		static string typeToCategory(const uint i);
		static bool isValidCategory(const string& c);
		static vector<uint> ensureCategoryToType(const vector<string>& c);
		static uint getNumCategories() { return _catmap.size(); }
		static void clearCategories();

	private:
		// labels -> 0,1,2,...,s
		// static: same across multiple sequences
		static map<string, uint> _catmap;
		static map<uint, string> _catrmap;
	};
}

#endif
