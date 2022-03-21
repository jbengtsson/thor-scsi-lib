#include <tps/utils.h>
#include <exception>
#include <iostream>
#include <limits.h>

double eps_tps  = 1e-25;

long int fact(long int n)
{
	if(n == 0 || n == 1){
		return 1;
	} else if (n > 1) {
		unsigned long result = 1;
		for(long int i = 2; i <= n; ++i){
			result *= i;
		}
		return result;
	} else {
		std::cout << "fact: neg. argument: " << n << "\n";
		throw std::domain_error("fact: neg. argument!");
	}
}

unsigned long binom(unsigned long n, unsigned long k)
{
	unsigned long c = 1, i;

	/*
	 * add checks if n < k ...
	 */
	if(n < k){
		return 0;
	}

	if (k > n-k){
		// take advantage of symmetry
		k = n-k;
	}

	for (i = 1; i <= k; i++, n--) {
		const unsigned long threshold = (n==0) ? ULONG_MAX : ULONG_MAX/n;
		if (c/i >= threshold){
			// raise an exception
			// assert(0);
			throw std::domain_error("binom c/l >= ULONG_MAX/n = threshold");
			return 0;
		}
		// split c * n / i into (c / i * i + c % i) * n / i
		c = c / i * n + c % i * n / i;
	}

	return c;
}


/**
 * @brief n over k
 *
 * @todo replace with
 */
long int nok_hist(long int n, long int k)
{
  long int j;
  double   u;

  u = 1.0;
  for (j = 0; j < k; j++)
    u *= (double)(n-j)/(double)(k-j);
  return (long int)(u+0.5);
}
