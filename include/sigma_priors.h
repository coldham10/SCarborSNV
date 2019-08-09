#ifndef SIGMA_PRIORS_H_
#define SIGMA_PRIORS_H_

typedef struct {
    long double l_P_tree; /*Calculated in sigma_priors.c */
    double c_thresh; /*Log of threshold for sigma=0 posterior to call candidates*/
    double l_lambda; /*Somatic mutation frequency*/
    double l_mu; /*Germline mutation frequency*/
    double l_P_H; /*Haploid event frequency. 0.09? this should be from empirical */
    /*FIXME get initial values for these two*/
    /*TODO add getopt options to change these following two*/
    double l_P_amp_err;
    double l_P_ADO;
    double l_P_clonal; /*0.51(?) refers to proportion of somatic genetic mutations that are public */
    int m;
} prior_params_t;

/* Returns the log prior probabilities for alternate allele count (sigma)
 * across a collection of m cells. Returned vector contains 2m+1 values
 * since sigma can range between 0 and 2m inclusive.
 */
int log_sigma_priors(prior_params_t* p, long double* l_priors);



#endif /*SIGMA_PRIORS_H_ */
