#include <stdlib.h>
#include "posteriors.h"
#include "math_utils.h"

int sigma_posteriors(long double* posteriors, long double* priors, long double* likelihoods, int m) {
    int i;
    long double sum;
    long double* products = malloc((2*m + 1) * sizeof(long double));
    for (i = 0; i < 2*m + 1; i++) {
        products[i] = priors[i] + likelihoods[i];
    }
    sum = LSE(products, 2*m + 1);
    for (i = 0; i < 2*m + 1; i++) {
        posteriors[i] = products[i] - sum;
    }
    free(products);
    return 0;
}

    
