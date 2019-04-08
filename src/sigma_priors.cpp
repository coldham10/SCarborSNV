#include <cstdlib>
#include <cmath>
#include "sigma_priors.h"

double log_prior(int m, int i);
long double log_factorial(int x);

std::vector<double> log_sigma_priors(int m) {
    std::vector<double> l_prior_vec; 
    for (int i = 0; i <= 2*m; i++) {
        l_prior_vec.push_back(log_prior(m, i));
    }
    return l_prior_vec;
}

double log_prior(int m, int i) {
    //FIXME
    return static_cast<double>(log_factorial(i));
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




    
