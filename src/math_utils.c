#include <math.h>
#include "math_utils.h"

/*Log-sum-exp */
long double LSE(long double* to_sum, int n) {
    int i;
    long double maxv = logl(0.0);
    long double sum = 0;
    for (i=0; i<n; i++) {
        maxv = (to_sum[i] > maxv) ? to_sum[i] : maxv;
    }
    if (maxv == logl(0.0)) {
        return maxv;
    }
    for (i=0; i<n; i++) {
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
