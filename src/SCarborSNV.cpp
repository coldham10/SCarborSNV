#include <iostream>
#include <cmath>
#include "sigma_priors.h"

int main(int argc, char** argv) {
    int m;
    if (argc < 2 ) {
        std::cout << "Enter m: ";
        std::cin >> m;
    }
    else {
        m = atoi(argv[1]);
    }
    double lambda0 = 0.0001;
//    double lambda0 = 0.01;
    double mu0 = 0.1;
    double P_H0 = 0.09;
    double P_clonal0 = 0.51;

    prior_params* p0 = new prior_params;
    p0->l_lambda    = std::log(lambda0);
    p0->l_mu        = std::log(mu0);
    p0->l_P_H       = std::log(P_H0);
    p0->l_P_clonal  = std::log(P_clonal0);
    p0->m           = m;

    double sum = 0;
    int i=0;
    std::vector<double> l_P_sigmas = log_sigma_priors(p0);
    for (std::vector<double>::iterator it = l_P_sigmas.begin(); it != l_P_sigmas.end(); ++it) {
        std::cout << i++ << ": " << std::exp(*it) << std::endl;
        sum += std::exp(*it);
    }
    std::cout <<"Sum: " << sum << std::endl;
    return 0;
}
