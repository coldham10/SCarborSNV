#include <stdlib.h>
#include <math.h>
#include "math_utils.h"

#define FREE_MEM (-123)

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
    else if (x == FREE_MEM) {
        /*signal to free memory */
        if (l_facts != NULL) {
            free(l_facts);
            l_facts = NULL;
        }
        return 0;
    }
    return l_facts[x];
}

long double l_binom_coeff(int n, int k) {
    return log_factorial(n) - log_factorial(k) - log_factorial(n-k);
}

long double l_binom_dist2(int m, int g, int sig) {
    static int last_m = -1;
    static long double* values = NULL;
    int s, j;
    long double log_p, log_q;
    if(m == FREE_MEM) {
        if (values == NULL) {
            /*data already freed or not yet alloc'd*/
            return 1;
        }
        else {
            last_m = -1;
            free(values);
            values = NULL;
            return 0;
        }
    }
    if(m != last_m) {
        /*Cannot use memoized values*/
        if (values != NULL) {
            free(values);
        }
        values = malloc((2*m + 1) * 3 * sizeof(long double));
        values[0] = logl(1);
        values[1] = values[2] = logl(0);
        for (s = 1; s < 2*m; s++) {
            log_p = logl(s) - logl(2*m);
            log_q = logl(1 - expl(log_p));
            for (j = 0; j < 3; j++) {
                values[3*s + j] = logl(1 + (j % 2)) + (j * log_p) + ((2-j) * log_q); 
            }
        }
        values[3*2*m + 0] = values[3*2*m + 1] = logl(0);
        values[3*2*m + 2] = logl(1);
        last_m = m;
    }
    return values[sig*3 + g];
}


void free_log_factorials() {
    /*Special signal*/
    log_factorial(FREE_MEM);
    return;
}

void free_log_binom() {
    l_binom_dist2(FREE_MEM, -1, -1);
    return;
}

