#include <iostream>
#include <cmath>
#include <getopt.h>
#include "scarborsnv.h"
#include "sigma_priors.h"

void init_params(global_params* gp, prior_params* p0, int argc, char** argv);

int main(int argc, char** argv) {
    //TODO eventually m will be calculated from number of BAMS/pileup info.
    int m;
    std::cout << "Enter m: ";
    std::cin >> m;
    prior_params p0{}; //*//= new prior_params;
    global_params gp{};//*// = new global_params;
    p0.m = m;
    init_params(&gp, &p0, argc, argv);

    double sum = 0;
    int i=0;
    std::vector<double> l_P_sigmas = log_sigma_priors(&p0);
    for (std::vector<double>::iterator it = l_P_sigmas.begin(); it != l_P_sigmas.end(); ++it) {
        std::cout << i++ << ": " << std::exp(*it) << std::endl;
        sum += std::exp(*it);
    }
    //delete p0;
    //delete gp;
    std::cout <<"Sum: " << sum << std::endl;
    return 0;
}

void init_params(global_params* gp, prior_params* p0, int argc, char** argv) {
    //Default values
    unsigned int n_threads = 1;
    std::string mp_fname = "";
    unsigned int mp_isfile = 0;
    double lambda0 = 0.0001; //0.01;
    double mu0 = 0.1;
    double P_H0 = 0.09;
    double P_clonal0 = 0.51;

    static struct option prior_options[] = {
        {"lambda", required_argument, NULL, 'l'},
        {"mu", required_argument, NULL, 'm'},
        {"p_haploid", required_argument, NULL, 'h'},
        {"p_clonal", required_argument, NULL, 'c'},
        {"pileup_file", required_argument, NULL, 'p'}
    };

    int c, opt_idx = 0;
    while ((c = getopt_long(argc, argv, "t:p:", prior_options, &opt_idx) )!= -1 ) {

        switch(c) {
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'p':
                mp_isfile = 1;
                mp_fname = optarg;
                break;
            case 'l':
                lambda0 = atof(optarg);
                break;
            case 'm':
                mu0 = atof(optarg);
                break;
            case 'h':
                P_H0 = atof(optarg);
                break;
            case 'c':
                P_clonal0 = atof(optarg);
                break;
        }
}
           
    gp->n_threads   = n_threads;
    gp->mp_isfile   = mp_isfile;
    //FIXME learn to use c++ strings...
    //std::strcpy(gp->mp_fname, mp_fname);
    gp->mp_fname     = mp_fname;

    p0->l_lambda    = std::log(lambda0);
    p0->l_mu        = std::log(mu0);
    p0->l_P_H       = std::log(P_H0);
    p0->l_P_clonal  = std::log(P_clonal0);
}

