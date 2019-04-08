#ifndef SIGMA_PRIORS_H_
#define SIGMA_PRIORS_H_

#include <vector>
typedef struct {
    double l_lambda;
    double l_P_H;
    unsigned int m;
} initial_consts;

/* Returns the log prior probabilities for alternate allele count (sigma)
 * across a collection of m cells. Returned vector contains 2m+1 values
 * since sigma can range between 0 and 2m inclusive.
 */
std::vector<double> log_sigma_priors(int m);



#endif //SIGMA_PRIORS_H_
