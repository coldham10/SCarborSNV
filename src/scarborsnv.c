#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "scarborsnv.h"
#include "math_utils.h"
#include "sigma_priors.h"
#include "pileup_reader.h"
#include "cell_likelihoods.h"
#include "locus_likelihoods.h"
#include "posteriors.h"

#define LOCUS_BATCH_SIZE (10)

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv);
void update_pairwise_p(Locus* loci_batch, int batch_size, prior_params_t* p0, long double* sig_priors, long double** numerators, int** denominators);

int main(int argc, char** argv) {
    int i, m, n_loci_read;
    long double* P_sigma;
    FILE* instream;

    prior_params_t* p0 = malloc(sizeof(prior_params_t));
    global_params_t* gp = malloc(sizeof(global_params_t));

    /*Uses getopt to parse cl parameters */
    init_params(gp, p0, argc, argv);

    m = gp->m;
    /* Open the mpileup stream*/
    if (gp->mp_isfile) {
        instream = fopen(gp->mp_fname, "r");
    }
    else {
        instream = stdin;
    }

    /* Compute 2m+1 priors and place in P_sigma */
    P_sigma = malloc((2*m + 1) * sizeof(long double));
    log_sigma_priors(p0, P_sigma);

    /*Two mxm matrices for computing pairwise p-bar*/
    long double** p_bar_numerators = malloc(m * sizeof(long double*));
    for (i = 0; i < m; i++) {
        p_bar_numerators[i] = malloc(m * sizeof(long double));
    }
    /*XXX remember to update both symmetic values */
    int** p_bar_denominators = malloc(m * sizeof(int*));
    for (i = 0; i < m; i++) {
        p_bar_denominators[i] = malloc(m * sizeof(int));
    }

    /* Retrieve & process loci in batches */
    Locus* loci_batch = malloc(LOCUS_BATCH_SIZE * sizeof(Locus));
    while (1) {
        n_loci_read = read_batch_loci(instream, loci_batch, LOCUS_BATCH_SIZE, m);
        if (n_loci_read == 0) {
            break;
        }
        update_pairwise_p(loci_batch, n_loci_read, p0, P_sigma, p_bar_numerators, p_bar_denominators);

        delete_locus_contents(loci_batch, n_loci_read, m);
    }

    if (gp->mp_isfile) {
        fclose(instream);
    }
    free(p0); free(gp);
    free(P_sigma);
    free(loci_batch);
    for (i = 0; i < m; i++) {
        free(p_bar_numerators[i]);
    }
    free(p_bar_numerators);
    for (i = 0; i < m; i++) {
        free(p_bar_denominators[i]);
    }
    free(p_bar_denominators);
    free_log_factorials();

    return 0;
}

void update_pairwise_p(Locus* loci_batch, int batch_size, prior_params_t* p, long double* sig_priors, long double** numerators, int** denominators) {
    int i, j;
    long double** cell_ls  = malloc(batch_size * sizeof(long double*));
    /*Posterior sigma distribution for each locus in batch*/
    long double** l_P_sig__D = malloc(batch_size * sizeof(long double*));
    long double*  locus_ls = malloc((2 * p->m + 1) * sizeof(long double));
    /*Iterate through loci for cell likelihoods and posterior sigma distributions*/
    for (i = 0; i < batch_size; i++) {
        cell_ls[i] = malloc(3 * p->m * sizeof(long double));
        for (j = 0; j < p->m; j++) {
            cell_likelihoods(&(loci_batch[i].cells[j]),
                    cell_ls[i] + 3*j,
                    loci_batch[i].ref_base,
                    p->l_P_amp_err,
                    p->l_P_ADO);
        }
        /*P(D|sigma) into locus_ls*/
        locus_likelihoods(cell_ls[i], locus_ls, p);
        /*P(sigma|D) into l_P_sig__D*/
        l_P_sig__D[i] = malloc((2*p->m + 1) * sizeof(long double));
        sigma_posteriors(l_P_sig__D[i], sig_priors, locus_ls, p->m);
    }
    free(locus_ls);
    for (i = 0; i < batch_size; i++) {
        for (j = 0; j < 2*p->m + 1; j++) {
            printf("%Lf ", expl(l_P_sig__D[i][j]));
        }
        printf("\n");
        free(l_P_sig__D[i]);
        free(cell_ls[i]);
    }
    free(l_P_sig__D);
    free(cell_ls);
    return;
}

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv) {
    /*Default values */
    unsigned int n_threads = 1;
    int m_cells = 1;
    int mp_isfile = 0;
    double lambda0 = 0.0001; 
    double mu0 = 0.1;
    double P_H0 = 0.09;
    double P_clonal0 = 0.51;
    /*These initial values are taken from the monovar source*/
    double P_amplification_err = 0.002;
    double P_ADO = 0.2;


    /*Using getopt to get command line arguments */
    /*TODO: handle getting bams */
    static struct option prior_options[] = {
        {"lambda", required_argument, NULL, 'l'},
        {"mu", required_argument, NULL, 'u'},
        {"p-haploid", required_argument, NULL, 'h'},
        {"p-clonal", required_argument, NULL, 'c'},
        {"pileup-file", required_argument, NULL, 'p'},
        {"n-cells", required_argument, NULL, 'm'}
    };
    int c, opt_idx = 0;
    /*TODO add options for amplification err and P_ADO*/
    while ((c = getopt_long(argc, argv, "t:p:m:", prior_options, &opt_idx) )!= -1 ) {
        switch(c) {
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'p':
                /*If pileup comes from file, -p pfile.pileup or --pileup-file=pfile.xx */
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
    p0->l_P_amp_err = log(P_amplification_err);
    p0->l_P_ADO     = log(P_ADO);
    p0->m           = m_cells;
}

