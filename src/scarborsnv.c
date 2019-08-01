#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "scarborsnv.h"
#include "sigma_priors.h"
#include "pileup_reader.h"

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv);

int main(int argc, char** argv) {
    int m;
    long double* P_sigma;
    prior_params_t* p0 = malloc(sizeof(prior_params_t));
    global_params_t* gp = malloc(sizeof(global_params_t));

    /*Uses getopt to parse cl parameters */
    init_params(gp, p0, argc, argv);

    m = gp->m;
    P_sigma = malloc((2*m + 1) * sizeof(long double));
    log_sigma_priors(p0, P_sigma);

    printf("priors computed\n");

    free(p0); free(gp);
    free(P_sigma);

    return 0;
}

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv) {
    /*Default values */
    unsigned int n_threads = 1;
    int m_cells = 1;
    /*char mp_fname[255] = ""; */
    int mp_isfile = 0;
    double lambda0 = 0.0001; /*0.01; */
    double mu0 = 0.1;
    double P_H0 = 0.09;
    double P_clonal0 = 0.51;

    /*Using getopt to get command line arguments */
    /*TODO: handle getting bams */
    static struct option prior_options[] = {
        {"lambda", required_argument, NULL, 'l'},
        {"mu", required_argument, NULL, 'u'},
        {"p_haploid", required_argument, NULL, 'h'},
        {"p_clonal", required_argument, NULL, 'c'},
        {"pileup-file", required_argument, NULL, 'p'},
        {"n-cells", required_argument, NULL, 'm'}
    };

    int c, opt_idx = 0;
    while ((c = getopt_long(argc, argv, "t:p:m:", prior_options, &opt_idx) )!= -1 ) {

        switch(c) {
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'p':
                /*If pileup comes from file, -p pfile.pileup or --pileup-file pfile.xx */
                mp_isfile = 1;
                strcpy(gp->mp_fname, optarg);
                break;
            case 'l':
                lambda0 = atof(optarg);
                break;
            case 'u':
                mu0 = atof(optarg);
                break;
            case 'h':
                P_H0 = atof(optarg);
                break;
            case 'c':
                P_clonal0 = atof(optarg);
                break;
            case 'm':
                m_cells = atoi(optarg);
                break;
        }
    }
           
    /*Populate global parameter struct */
    gp->n_threads   = n_threads;
    gp->mp_isfile   = mp_isfile;
    gp->m           = m_cells;

    /*Convert prior params to log space, populate intial parameters struct */
    p0->l_lambda    = log(lambda0);
    p0->l_mu        = log(mu0);
    p0->l_P_H       = log(P_H0);
    p0->l_P_clonal  = log(P_clonal0);
    p0->m           = m_cells;
}

