#ifndef __UTILS_H__
#define __UTILS_H__

#include "params.hpp"

#include <map>
#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <utility>
#include <sstream>
#include <iomanip>
#include <list>
using namespace std;

namespace HMM_PP {

	typedef pair<ullint, ullint> ullintPair;

	/*
	 * Reads non-white-space lines from stream and inputs them to strings.
	 * Note that stream should NOT be accessed externally once passed to istreamLineReader object.
	 */
	class istreamLineReader {
	public:
		istreamLineReader(istream& stream, bool deleteInner = false);
		~istreamLineReader();

		bool eof();
		istreamLineReader& operator>>(string& s);

	private:
		istream& _innerStream;
		bool _deleteInner;

		string _nextToExtract;

		void advanceToNonBlankLine();

		static const string BLANK;
	};

	class ostreamWriter {
	public:
		ostreamWriter(ostream& stream, bool deleteInner = false);
		~ostreamWriter();

		inline ostream& operator()() { return _stream; }

	private:
		ostream& _stream;
		bool _deleteInner;
	};

	class utils {
	public:
		static const string STD_STREAM;

		static istreamLineReader* getIstreamLineReaderFromFile(string file, bool requireNamedFile = false);
		static ostreamWriter* getOstreamWriterFromFile(string file);
	};

	template<class T>
	void size(vector<vector<T> >& d, const uint r, const uint c) {
		d.clear();
		d.resize(r);
		for (uint i = 0; i < r; i++)
			d[i].resize(c);
	}

	template<class T>
	inline T sign(T t) {
		return (t < 0) ? T(-1) : T(1);
	}

	template<class T>
	inline T& capAbsoluteValue(T& t, const T& capVal) {
		if (abs(t) > capVal)
			t = sign(t) * capVal;

		return t;
	}

	template<class T>
	class PrecisionThreshold {
	public:
		PrecisionThreshold(int precision = 0, T thresh = 0) : _precision(precision), _thresh(thresh) {}
		~PrecisionThreshold() {}

		bool passScore(T eval) const {
			setPrecision(eval);
			return eval >= _thresh;
		}

		bool failScore(T eval) const {
			setPrecision(eval);
			return eval < _thresh;
		}

		friend ostream& operator<<(ostream& stream, const PrecisionThreshold& t) {
			T tmpThresh = t._thresh;
			t.setPrecision(tmpThresh);
			stringstream str;
			str << tmpThresh;
			stream << str.str();

			return stream;
		}

		const T& getThreshold() const { return _thresh; }

	private:
		int _precision;
		T _thresh;

		void setPrecision(T& eval) const {
			stringstream stream;
			stream << setiosflags(ios::fixed) << setprecision(_precision) << eval;
			stream >> eval;
		}
	};

	string int2str(const uint i);

	bool fileExists(const string& f);
	bool removeFile(const string& f);

	bool beginsWith(string str, string prefix);
	bool endsWith(string str, string suffix);

	// Math functions
	bool realNum(double);

	/*
	// returns log (e^logprob1 + e^logprob2).
	double sumLogProb(double logprob1, double logprob2);

	// The input array contains a set of log probabilities lp1, lp2, lp3
	// The return value should be the log of the sum of the
	//   probabilities: log(e^lp1 + e^lp2 + e^lp3 + ...)
	double sumLogProb(const vector<double>& logprobs);
	 */

	/*
	bool constrain(double*, double, double);
	 */

	template<typename T, size_t N>
	list<T> arrayToList(const T (&data)[N]) {
	    return list<T>(data, data+N);
	}
}

#endif
