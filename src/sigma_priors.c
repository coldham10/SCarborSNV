#include <stdlib.h>
#include <math.h>
#include "sigma_priors.h"
#include "math_utils.h"

/* Note on function names, for example: "l_P_sig__sSNV_NSNVT_H_NHT" is read as
 * log probability of sigma given there is a somatic SNV shared by all cells
 * and there is a haploid event shared by all cells */


/*************************************
 *
 * Counting Tools
 *
 * *********************************/

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
    else if (x == -123) {
        /*signal to free memory */
        free(l_facts);
        l_facts = NULL;
        return -123.0;
    }
    return l_facts[x];
}

long double log_binom(int n, int k) {
    return log_factorial(n) - log_factorial(k) - log_factorial(n-k);

}


/*************************************
 *
 * Mutant Priors
 *
 * *********************************/


/*"T" function in paper */
long double l_T(int a, int b) {
    /*Kuipers et al. inspired tree function */
    long double numerator = 2 * log_binom(a,b);
    long double denominator = logl(2*b-1) + log_binom(2*a, 2*b);
    return numerator - denominator;
}
    

long double l_P_ancestral__subclonal(int m) {
    /*The log probability that a given subclonal mutation is ancestral to all sequenced cells. */
    /*Assumes a neutral evolution model. */
    if (m >= 50) {
        return logl(0.0);
    }
    long double f_lower = 1.0e-8;
    long double f_upper = 0.5;
    /*By integrating P(f): */
    long double l_k = -1.0 * logl((logl(f_upper/f_lower) - 2*f_upper + 2*f_lower));

    /*TODO check for underflow and overflow */
    /*log[upper^m - lower^m] */
    long double diff_pow_m  = logl(powl(f_upper, m) - powl(f_lower, m));
    /*log[upper^(m+1) - lower^(m+1)] */
    long double diff_pow_m1  = logl(powl(f_upper, m+1) - powl(f_lower, m+1));
    /*first term in integration */
    long double term_1 = m*logl(2.0) - logl(m);
    term_1 = term_1 + diff_pow_m;
    /*second term in integration */
    long double term_2 = (m+1)*logl(2.0) - logl(m+1);
    term_2 = term_2 + diff_pow_m1;

    return l_k + logl(expl(term_1) - expl(term_2));
}

long double l_P_intree__somatic(int m, long double l_P_clonal) {
    long double P_clonal = expl(l_P_clonal);
    long double l_P_ancestral = l_P_ancestral__subclonal(m) + logl(1-P_clonal);
    return logl(1 - P_clonal - expl(l_P_ancestral));
}

long double l_P_sig__sSNV_NSNVT_H_NHT(int m, int sig) {
     return logl((long double)((sig == 0 || sig == 2*m) ? 0.5 : 0 ));
}

long double l_P_sig__sSNV_NSNVT_HT(int m, int sig) {
    if (sig > m) {
        return logl(0.5) + l_T(m, sig-m);
    } else if (m > sig) {
        return logl(0.5) + l_T(m, m-sig);
    } else {
        return logl(0.0);
    }
}

long double l_P_sig__sSNV_NSNVT_H(int m, int sig, long double l_P_HT__H) {
    long double term_1 = logl(1-expl(l_P_HT__H)) + l_P_sig__sSNV_NSNVT_H_NHT(m, sig);
    long double term_2 = l_P_HT__H + l_P_sig__sSNV_NSNVT_HT(m, sig);
    return LSE2(term_1, term_2);
}

long double l_P_sig__sSNV_NSNVT_NH(int m, int sig) {
     return logl((long double)((sig == m) ? 1 : 0 ));
}

long double l_P_sig__sSNV_NSNVT(int m, int sig, long double l_P_H, long double l_P_HT__H) {
    long double term_1 = l_P_sig__sSNV_NSNVT_H(m, sig, l_P_HT__H) + l_P_H;
    long double term_2 = l_P_sig__sSNV_NSNVT_NH(m, sig) + logl(1-expl(l_P_H));
    return LSE2(term_1, term_2);
}

long double l_P_sig__SNVT_HT(int m, int sig) {
    int a, h;
    int n_sum = 0;
    int array_length = (int)(m/2 + 1);
    long double* to_sum = malloc(array_length*sizeof(long double));
    for (a = 1; a < m; ++a) {
        for (h = 1; h < m; ++h) {
            if (n_sum + 3 > array_length) {
                array_length *= 2;
                to_sum = realloc(to_sum, array_length*sizeof(long double));
            }
            if (a-h == sig && h <= a) {
                /*Case 1 */
                to_sum[n_sum++] = l_T(m,a) + l_T(a,h);
            }
            if (a+h == sig && h <= a) {
                /*Case 2 */
                to_sum[n_sum++] = l_T(m,a) + l_T(a,h);
            }
            if (2*a == sig && h >= a) {
                /*Case 3 */
                to_sum[n_sum++] = l_T(m,h) + l_T(h,a);
            }
        } /*for h */
    } /*for a */
    long double l_sum = LSE(to_sum, n_sum);
    free(to_sum);
    return l_sum - logl(3) + logl(2*m-1) - logl(2*(m-1));
}

long double l_P_sig__SNVT_H_NHT(int m, int sig) {
    if (sig > 0 && sig < 2*m && sig % 2 == 0) {
        return logl(2*m - 1) - logl(2*(m-1)) + l_T(m, sig/2);
    }
    else {
        return logl(0);
    }
}

long double l_P_sig__SNVT_H(int m, int sig, long double l_P_HT__H) {
    long double term_1 = l_P_sig__SNVT_HT(m, sig) + l_P_HT__H;
    long double term_2 = l_P_sig__SNVT_H_NHT(m, sig) + logl(1-expl(l_P_HT__H));
    return LSE2(term_1, term_2);
}

long double l_P_sig__SNVT_NH(int m, int sig) {
    if (sig == 0 || sig >= m) {
        return logl(0);
    }
    return logl(2*m - 1) - logl(2*(m-1)) + l_T(m, sig);
}

long double l_P_sig__SNVT(int m, int sig, long double l_P_H, long double l_P_HT__H) {
    long double term_1 = l_P_sig__SNVT_H(m, sig, l_P_HT__H) + l_P_H;
    long double term_2 = l_P_sig__SNVT_NH(m, sig) + logl(1-expl(l_P_H));
    return LSE2(term_1, term_2);
}

long double l_P_sig__sSNV(int m, int sig, long double l_P_H, long double l_P_tree) {
    /*FIXME does not sum to 1 if m=1 */
    long double term_1 = l_P_sig__SNVT(m, sig, l_P_H, l_P_tree) + l_P_tree;
    long double term_2 = l_P_sig__sSNV_NSNVT(m, sig, l_P_H, l_P_tree) + logl(1-expl(l_P_tree));
    return LSE2(term_1, term_2);
}

/*************************************
 *
 * Welltype Priors
 *
 * *********************************/


long double l_P_sig__NsSNV_H(int m, int sig, long double l_mu, long double l_P_HT__H) {
    long double term_1, term_2;
    long double l_om_mu = logl(1-expl(l_mu));
    long double l_om_P_HT__H = logl(1-expl(l_P_HT__H));
    if (sig == 0) {
        term_1 = 2 * l_om_mu;
        term_2 = l_mu + l_om_mu + l_om_P_HT__H;
        return LSE2(term_1, term_2);
    } else if (sig > 0 && sig < 2*m && sig != m) {
        term_1 = l_mu + l_om_mu + l_P_HT__H;
        term_2 = l_T(m, abs(sig-m)) + logl(2*m-1) - logl(2*(m-1));
        return term_1 + term_2;
    } else if (sig == m) {
        return logl(0);
    } else if (sig == 2*m) {
        term_1 = 2 * l_mu;
        term_2 = l_mu + l_om_mu + l_om_P_HT__H;
        return LSE2(term_1, term_2);
    } else {
        return -1;
    }
}

long double l_P_sig__NsSNV_NH(int m, int sig, long double l_mu) {
    /*HWE */
    long double l_om_mu = logl(1-expl(l_mu));
    if (sig == 0) {
        return 2*l_om_mu;
    } else if (sig == m) {
        return logl(2) + l_mu + l_om_mu;
    } else if (sig == 2*m) {
        return 2*l_mu;
    } else {
        return logl(0);
    }
}

long double l_P_sig__NsSNV(int m, int sig, long double l_mu, long double l_P_H, long double l_P_HT__H) {
    long double term_1 = l_P_sig__NsSNV_H(m, sig, l_mu, l_P_HT__H) + l_P_H;
    long double term_2 = l_P_sig__NsSNV_NH(m, sig, l_mu) + logl(1-expl(l_P_H));
    return LSE2(term_1, term_2);
}

long double l_P_sig(int m, int sig, prior_params_t* p) {
    long double l_om_lambda = logl(1-expl(p->l_lambda));
    long double term_1 = p->l_lambda + l_P_sig__sSNV(m, sig, p->l_P_H, p->l_P_tree);
    long double term_2 = l_om_lambda + l_P_sig__NsSNV(m, sig, p->l_mu, p->l_P_H, p->l_P_tree);
    return LSE2(term_1, term_2);
}

/*Function to get overall priors */
int log_sigma_priors(prior_params_t* p, long double* l_priors) {
    /*FIXME: does not sum to 1 if m=1 */
    int sig;
    int m = p->m;
    p->l_P_tree = l_P_intree__somatic(m, p->l_P_clonal);
    for (sig = 0; sig <= 2*m; sig++) {
        l_priors[sig] = l_P_sig(m, sig, p);

    }
    /*Clear factorial memory */
    log_factorial(-123);
    return 0;
}
