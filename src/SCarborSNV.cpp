#include <iostream>
#include <cmath>
#include "sigma_priors.h"

int main(void) {
    std::cout << "Enter m: ";
    int m;
    std::cin >> m;
    std::vector<double> l_P_sigmas = log_sigma_priors(m);
    for (std::vector<double>::iterator it = l_P_sigmas.begin(); it != l_P_sigmas.end(); ++it) {
        std::cout << std::exp(*it) << std::endl;
    }
    return 0;
}
