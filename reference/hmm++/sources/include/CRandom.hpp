#ifndef __CRANDOM_H__
#define __CRANDOM_H__

#include <vector>
using namespace std;

class CRandom {
public:

	static const int IA;
	static const int IM;
	static const int IQ;
	static const int IR;
	static const int NTAB;
	static const int NDIV;

	static const double EPS;
	static const double AM;
	static const double RNMX;

	// Current seed
	static int idum;

	static int iy;
	static vector<int> iv;

	static double last;

	static void srand(long unsigned iseed = 0);
	static double rand();
	static int rand(int);

};

#endif
