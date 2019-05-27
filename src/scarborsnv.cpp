#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include "scarborsnv.h"
#include "piler.h"
#include "sigma_priors.h"

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv);

int main(int argc, char** argv) {
    prior_params_t p0{};
    global_params_t gp{};
    init_params(&gp, &p0, argc, argv);

    piler_module::Piler* piler;
    if (gp.mp_isfile) {
        std::ifstream ifs(gp.mp_fname, std::ifstream::in);
        piler = new piler_module::Piler(&ifs, gp.n_threads);
    }
    else {
        piler = new piler_module::Piler(&std::cin, gp.n_threads);
    }


    delete piler;
    return 0;
}

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv) {
    //Default values
    unsigned int n_threads = 1;
    std::string mp_fname = "";
    bool mp_isfile = false;
    double lambda0 = 0.0001; //0.01;
    double mu0 = 0.1;
    double P_H0 = 0.09;
    double P_clonal0 = 0.51;

    //Using getopt to get command line arguments
    //TODO: handle getting bams
    static struct option prior_options[] = {
        {"lambda", required_argument, NULL, 'l'},
        {"mu", required_argument, NULL, 'm'},
        {"p_haploid", required_argument, NULL, 'h'},
        {"p_clonal", required_argument, NULL, 'c'},
        {"pileup-file", required_argument, NULL, 'p'}
    };

    int c, opt_idx = 0;
    while ((c = getopt_long(argc, argv, "t:p:", prior_options, &opt_idx) )!= -1 ) {

        switch(c) {
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'p':
                //If pileup comes from file, -p pfile.pileup or --pileup-file pfile.xx
                mp_isfile = true;
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
           
    //Populate global parameter struct
    gp->n_threads   = n_threads;
    gp->mp_isfile   = mp_isfile;
    gp->mp_fname    = mp_fname;

    //Convert prior params to log space, populate intial parameters struct
    p0->l_lambda    = std::log(lambda0);
    p0->l_mu        = std::log(mu0);
    p0->l_P_H       = std::log(P_H0);
    p0->l_P_clonal  = std::log(P_clonal0);
    //FIXME a hack to insert m
    p0->m = n_threads;
}

