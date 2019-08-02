#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

/*Log-sum-exp using the 'trick' to prevent underflow*/
long double LSE(long double* to_sum, int n);
/*Two argument case of the above*/
long double LSE2(long double arg1, long double arg2);

#endif /*MATH_UTILS_H_*/
