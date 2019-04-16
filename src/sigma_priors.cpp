#include <cstdlib>
#include <cmath>
#include "sigma_priors.h"

double log_prior(int m, int sig);
long double log_factorial(int x);
long double l_P_intree__somatic(int m, double l_P_clonal);
//long double l_P_sig__sSNV_NSNVT(int m, int sig, double l_P_H, long double l_P_HT__H);
long double l_P_sig__SNVT_H(int m, int sig, long double l_P_HT__H);

std::vector<double> log_sigma_priors(int m) {
    std::vector<double> l_prior_vec; 
    //todo input prior params struct

    for (int i = 0; i <= 2*m; i++) {
        l_prior_vec.push_back(log_prior(m, i));
    }
    return l_prior_vec;
}

double log_prior(int m, int sig) {
    //FIXME
    static long double l_P_tree = l_P_intree__somatic(m, std::log(0.51));
    //return l_P_sig__sSNV_NSNVT(m, sig, std::log(0.09), l_P_tree);
    return l_P_sig__SNVT_H(m, sig, l_P_tree);
}

long double log_factorial(int x) {
    static int biggest = 0;
    static std::vector<long double> l_facts{static_cast<long double>(0)}; //ln(0!) = 0
    if (x > biggest) {
        for (int i = biggest + 1; i <= x; i++) {
            l_facts.push_back(std::log(static_cast<long double>(i)) + l_facts[i-1]);
        }
        biggest = x;
    }
    return l_facts[x];
}

inline long double log_binom(int n, int k) {
    return log_factorial(n) - log_factorial(k) - log_factorial(n-k);

}

long double LSE(std::vector<long double> to_sum) {
    long double maxv = std::log(static_cast<long double>(0));
    long double sum = 0;
    for (std::vector<long double>::iterator it = to_sum.begin(); it != to_sum.end(); ++it) {
        maxv = (*it > maxv) ? *it : maxv;
    }
    if (maxv == std::log(static_cast<long double>(0))) {
        return maxv;
    }
    for (std::vector<long double>::iterator it = to_sum.begin(); it != to_sum.end(); ++it) {
        //LSE 'trick'
        sum += std::exp(*it - maxv);
    }
    return std::log(sum) + maxv;
}

inline long double LSE(long double arg1, long double arg2) {
    std::vector<long double> v{arg1, arg2};
    return LSE(v);
}
long double l_T(int a, int b) {
    //Kuipers et al. inspired tree function
    long double numerator = 2 * log_binom(a,b);
    long double denominator = std::log(static_cast<long double>(2*b-1)) + log_binom(2*a, 2*b);
    return numerator - denominator;
}
    
long double l_P_ancestral__subclonal(int m) {
    //The log probability that a given subclonal mutation is ancestral to all sequenced cells.
    //Assumes a neutral evolution model.
    //TODO check that 50 is sufficient
    if (m >= 50) {
        return std::log(static_cast<long double>(0));
    }
    long double f_lower = 1.0e-8;
    long double f_upper = 0.5;
    //By integrating P(f):
    long double l_k = -1.0 * std::log((std::log(f_upper/f_lower) - 2*f_upper + 2*f_lower));

    //TODO check for underflow and overflow
    //log[upper^m - lower^m]
    long double diff_pow_m  = std::log(std::pow(f_upper, m) - std::pow(f_lower, m));
    //log[upper^(m+1) - lower^(m+1)]
    long double diff_pow_m1  = std::log(std::pow(f_upper, m+1) - std::pow(f_lower, m+1));
    //first term in integration
    long double term_1 = (m)*std::log(static_cast<long double>(2)) - std::log(static_cast<long double>(m));
    term_1 = term_1 + diff_pow_m;
    //second term in integration
    long double term_2 = (m+1)*std::log(static_cast<long double>(2)) - std::log(static_cast<long double>(m+1));
    term_2 = term_2 + diff_pow_m1;

    //TODO inplement more underflow resistant function to handle this sort of situation.
    return l_k + std::log(std::exp(term_1) - std::exp(term_2));
}

long double l_P_intree__somatic(int m, double l_P_clonal) {
    long double P_clonal = std::exp(static_cast<long double>(l_P_clonal));
    long double l_P_ancestral = l_P_ancestral__subclonal(m) + std::log(1-P_clonal);
    //TODO logaddexp
    return std::log(1 - P_clonal - std::exp(l_P_ancestral));
}

inline long double l_P_sig__sSNV_NSNVT_H_NHT(int m, int sig) {
     return std::log(static_cast<long double>((sig == 0 || sig == 2*m) ? 0.5 : 0 ));
}

long double l_P_sig__sSNV_NSNVT_HT(int m, int sig) {
    if (sig > m) {
        return std::log(0.5) + l_T(m, sig-m);
    } else if (m > sig) {
        return std::log(0.5) + l_T(m, m-sig);
    } else {
        return std::log(0.0);
    }
}

long double l_P_sig__sSNV_NSNVT_H(int m, int sig, long double l_P_HT__H) {
    long double term_1 = std::log(1-std::exp(l_P_HT__H)) + l_P_sig__sSNV_NSNVT_H_NHT(m, sig);
    long double term_2 = l_P_HT__H + l_P_sig__sSNV_NSNVT_HT(m, sig);
    return LSE(term_1, term_2);
}

inline long double l_P_sig__sSNV_NSNVT_NH(int m, int sig) {
     return std::log(static_cast<long double>((sig == m) ? 1 : 0 ));
}

long double l_P_sig__sSNV_NSNVT(int m, int sig, double l_P_H, long double l_P_HT__H) {
    long double term_1 = l_P_sig__sSNV_NSNVT_H(m, sig, l_P_HT__H) + static_cast<long double>(l_P_H);
    long double term_2 = l_P_sig__sSNV_NSNVT_NH(m, sig) + std::log(1-std::exp(static_cast<long double>(l_P_H)));
    return LSE(term_1, term_2);
}



long double l_P_sig__SNVT_HT(int m, int sig) {
    std::vector<long double> to_sum;
    for (int a = 1; a < m; ++a) {
        for (int h = 1; h < m; ++h) {
            if (a-h == sig && h <= a) {
                //Case 1
                to_sum.push_back(l_T(m,a) + l_T(a,h));
            }
            if (a+h == sig && h <= a) {
                //Case 2
                to_sum.push_back(l_T(m,a) + l_T(a,h));
            }
            if (2*a == sig && h >= a) {
                //Case 3
                to_sum.push_back(l_T(m,h) + l_T(h,a));
            }
        } //for h
    } //for a
    long double l_sum = LSE(to_sum);
    return l_sum - std::log(static_cast<long double>(3)) + std::log(2*m-1) - std::log(2*(m-1));
}

long double l_P_sig__SNVT_H_NHT(int m, int sig) {
    if (sig > 0 && sig < 2*m && sig % 2 == 0) {
        return std::log(static_cast<long double>(2*m - 1)) 
            - std::log(static_cast<long double>(2*(m-1)))
            + l_T(m, sig/2);
    } else {
        return std::log(static_cast<long double>(0));
    }
}

long double l_P_sig__SNVT_H(int m, int sig, long double l_P_HT__H) {
    long double term_1 = l_P_sig__SNVT_HT(m, sig) + l_P_HT__H;
    long double term_2 = l_P_sig__SNVT_H_NHT(m, sig) + std::log(1-std::exp(l_P_HT__H));
    return LSE(term_1, term_2);
}
