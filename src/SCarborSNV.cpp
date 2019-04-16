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
    double sum = 0;
    std::vector<double> l_P_sigmas = log_sigma_priors(m);
    for (std::vector<double>::iterator it = l_P_sigmas.begin(); it != l_P_sigmas.end(); ++it) {
//        std::cout << std::exp(*it) << std::endl;
        sum += std::exp(*it);
    }
    std::cout << sum << std::endl;
    return 0;
}
