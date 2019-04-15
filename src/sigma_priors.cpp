#include <cstdlib>
#include <cmath>
#include "sigma_priors.h"

double log_prior(int m, int sig);
long double log_factorial(int x);
long double l_P_intree__somatic(int m, double l_P_clonal);

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
    return static_cast<double>(l_P_tree);
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
        return static_cast<long double>(0);
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
    //FIXME currently term 2 is larger than term 1, giving a negative value to log!
    return l_k + std::log(std::exp(term_1) - std::exp(term_2));
}

long double l_P_intree__somatic(int m, double l_P_clonal) {
    long double P_clonal = std::exp(static_cast<long double>(l_P_clonal));
    long double l_P_ancestral = l_P_ancestral__subclonal(m) + std::log(1-P_clonal);
    //TODO logaddexp
    return std::log(P_clonal + std::exp(l_P_ancestral));
}


