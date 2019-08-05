#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

/*Log-sum-exp using the 'trick' to prevent underflow*/
long double LSE(long double* to_sum, int n);
/*Two argument case of the above*/
long double LSE2(long double arg1, long double arg2);
/*Returns the log of binomial coefficient (n k)*/
long double l_binom_coeff(int n, int k);
/* Returns log probability of g from binomial distribution B(2, sig/2m)
 * Keeps static memory, call free_log_binom() */
long double l_binom_dist2(int m, int g, int sig);
/*Factorial values are memoized, so this function frees that memory*/
void free_log_factorials();
void free_log_binom();

#endif /*MATH_UTILS_H_*/
