#include <stdlib.h>
#include <math.h>
#include "math_utils.h"

#define FREE_FACTORIAL (-123)

/*Log-sum-exp */
long double LSE(long double* to_sum, int n) {
    int i;
    long double maxv = logl(0.0);
    long double sum = 0;
    for (i = 0; i < n; i++) {
        maxv = (to_sum[i] > maxv) ? to_sum[i] : maxv;
    }
    if (maxv == logl(0.0)) {
        return maxv;
    }
    for (i = 0; i < n; i++) {
        /*LSE 'trick' */
        sum += expl(to_sum[i] - maxv);
    }
    return logl(sum) + maxv;
}

/*Two argument case */
long double LSE2(long double arg1, long double arg2) {
    long double big = (arg1 > arg2) ? arg1 : arg2;
    long double sml = (arg1 > arg2) ? arg2 : arg1;
    if (sml == logl(0)) {
        return big;
    }
    else {
        return big + logl(1 + expl(sml-big));
    }
}

long double log_factorial(int x) {
    static int biggest = 0;
    static long double* l_facts = NULL;

    int i;
    if (l_facts == NULL) {
        /*initialize static array */
        l_facts = calloc(1, sizeof(long double)); /*ln(0!) = 0 */
    }
    if (x > biggest) {
        l_facts = realloc(l_facts, (x+8)*sizeof(long double));
        for (i = biggest + 1; i < x+8; i++) {
            l_facts[i] = l_facts[i-1] + logl((long double)i);
        }
        biggest = x+7;
    }
    else if (x == FREE_FACTORIAL) {
        /*signal to free memory */
        if (l_facts != NULL) {
            free(l_facts);
            l_facts = NULL;
        }
        return 0;
    }
    return l_facts[x];
}

long double log_binom(int n, int k) {
    return log_factorial(n) - log_factorial(k) - log_factorial(n-k);
}

void free_log_factorials() {
    /*Special signal*/
    log_factorial(FREE_FACTORIAL);
    return;
}
