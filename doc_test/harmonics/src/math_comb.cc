#include "math_comb.h"
#include <limits.h>
#include <cassert>

unsigned long thor_scsi::binom(unsigned long n, unsigned long k) {
	unsigned long c = 1, i;

	if (k > n-k) // take advantage of symmetry
		k = n-k;

	for (i = 1; i <= k; i++, n--) {
		if (c/i >= ULONG_MAX/n){
			// raise an exception
			assert(0);
			return 0;
		}
		// split c * n / i into (c / i * i + c % i) * n / i
		c = c / i * n + c % i * n / i;
	}

	return c;
}

#if 0
double thor_scsi::binom (int const n, int const k)
{
	int i;
	double num, den;

	if (k<0  || n<0 || k>n){
		return -1.0;
	}

	if (k==0 || k==n){
		return 1.0;
	}

	if (k==1 || k==n-1){
		return n;
	}

	num = den = 1.0;

	for(i=n-k+1; i<=n; i++){
		num *= i;
	}
	for(i=1;     i<=k; i++){
		den *= i;
	}

	return num/den;
}
#endif
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
