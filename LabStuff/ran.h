//#include "nr3.h" 
#include <math.h>
#ifndef PI
#define PI 3.14159265358979323846264338328
#endif

typedef unsigned long long int Ullong;
typedef double Doub;
typedef unsigned int Uint;

struct Ran {
/*Implementation of the highest quality recommended generator. The constructor is called with
an integer seed and creates an instance of the generator. The member functions int64, doub,
and int32 return the next values in the random sequence, as a variable type indicated by their
names. The period of the generator is 3:138 1057 . */

	Ullong u,v,w;
	double phi, uu, rr, xx, yy;
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
	//Constructor. Call with any integer seed (except value of v above).
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	
	//Return 64-bit random integer. See text for explanation of method.
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	
	//Return random double-precision floating value in the range 0. to 1.
	inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
	
	//Return 32-bit random integer.
	inline Uint int32() { return (Uint)int64(); }

	//RETURN with respect to the gaussian distribution
	//check again the book of Werner Krauth (SMAC)
	inline double gauss(double sigma){
		phi = 2.0*PI*doub();
		uu = -log(doub());
		rr = sigma*sqrt(2.0*uu);
		xx = rr*cos(phi);
		yy = rr*sin(phi);
		return xx;
	}
};

