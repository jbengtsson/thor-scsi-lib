#ifndef _TPS_UTILS_H_
#define _TPS_UTILS_H_

/**
 *
 * Todo:
 *     Replace with inline function,
 *     check if GSL' fcmp is a better solution
 *
 *
 * Warning:
 *     Does not return 0 if value is 0
 *     Not side effect free
 */
#define sgn(n)      ((n > 0) ? 1 : ((n < 0) ? -1 : 0))

long int nok_using_dbl(long int n, long int k);

/**
 * @brief calculate binomial coefficient
 *
 * @ todo: check if in standard library or in GSL
 */
unsigned long binom(unsigned long n, unsigned long k);

inline long int nok(long int n, long int k){
	return binom(n, k);
}

long int nok_hist(long int n, long int k);


long int fact(long int n);

#endif /* _TPS_UTILS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
