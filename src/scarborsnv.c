#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "scarborsnv.h"
#include "math_utils.h"
#include "sigma_priors.h"
#include "pileup_reader.h"
#include "likelihoods.h"
#include "posteriors.h"
#include "tree.h"
#include "inference.h"

#define LOCUS_BATCH_SIZE (10)

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv);
long double read_candidate(FILE* c_file, Candidate cand);
int update_candidates(Locus* loci_batch,
        int batch_size,
        prior_params_t* p0,
        long double* sig_priors,
        FILE* candidates,
        long double** numr,
        int** denom);

int main(int argc, char** argv) {
    int i, j, m, n_loci_read;
    long double* P_sigma;
    FILE* instream;
    FILE* f_candidates;
    long double** p_bar_numerators;
    long double** distance_matrix;
    int** p_bar_denominators;
    unsigned long candidates_found = 0;
    Node* T;
    Candidate* candidate;
    prior_params_t* p0 = malloc(sizeof(prior_params_t));
    global_params_t* gp = malloc(sizeof(global_params_t));
    /*Uses getopt to parse cl parameters */
    init_params(gp, p0, argc, argv);
    m = gp->m;
    /*Open temp file for storing candidate genotype posteriors*/
    f_candidates = fopen(gp->tmp_fname,"wb+");
    if (f_candidates == NULL) {
        free(p0); free(gp);
        fprintf(stderr, "Could not open temp file\n");
        abort();
    }
    /* Open the mpileup stream*/
    if (gp->mp_isfile) {
        instream = fopen(gp->mp_fname, "r");
        if (instream == NULL) {
            free(p0); free(gp);
            fclose(f_candidates);
            fprintf(stderr, "Could not open pileup file\n");
            abort();
        }
    }
    else {
        instream = stdin;
    }
    /* Compute 2m+1 sigma priors and place in P_sigma */
    P_sigma = malloc((2*m + 1) * sizeof(long double));
    log_sigma_priors(p0, P_sigma);
    /*Two (m+1)x(m+1) matrices for computing pairwise p-bar, (incl. root)*/
    p_bar_numerators = malloc((m + 1) * sizeof(long double*));
    for (i = 0; i < m + 1; i++) {
        p_bar_numerators[i] = calloc(m + 1, sizeof(long double));
        for (j = 0; j < m + 1; j++) {
            p_bar_numerators[i][j] = logl(0);
        }
    }
    p_bar_denominators = malloc((m + 1) * sizeof(int*));
    for (i = 0; i < m + 1; i++) { p_bar_denominators[i] = calloc(m + 1, sizeof(int)); }
    /* Retrieve & process loci in batches, identify and store candidates*/
    Locus* loci_batch = malloc(LOCUS_BATCH_SIZE * sizeof(Locus));
    while (1) {
        n_loci_read = read_batch_loci(instream, loci_batch, LOCUS_BATCH_SIZE, m);
        if (n_loci_read == 0) {
            break;
        }
        /*Update p_bar matrices and add candidates to temp file*/
        candidates_found += update_candidates(loci_batch,
                n_loci_read,
                p0,
                P_sigma,
                f_candidates,
                p_bar_numerators,
                p_bar_denominators);
        delete_locus_contents(loci_batch, n_loci_read, m);
    }
    free(loci_batch);
    /*Done reading pileup file*/
    if (gp->mp_isfile) { fclose(instream); }
    fprintf(stderr, "Found %ld candidate loci\n", candidates_found);
    fprintf(stderr, "%d pairwise comparisons made between cells to infer phylogeny\n", sqr_mat_sum(p_bar_denominators, m)/2);
    /*Compute tree distances*/
    distance_matrix = malloc((m + 1) * sizeof(long double*));
    for (i = 0; i < m + 1; i++) { distance_matrix[i] = malloc((m + 1) * sizeof(long double)); }
    /*TODO orphan unsupported cells, update raw dists, m, etc*/
    i = expected_jukes_cantor(distance_matrix, p_bar_numerators, p_bar_denominators, m + 1);
    if (i != 0) {
        fprintf(stderr, "Failed to find all pairwise distances. Some cells may share no overlapping loci\n");
         /* If still nans after orphaning e.g. "Fast NJ-like algorithms to deal with incomplete distance matrices" --google.*/
    }
    /*Freeing old matrices*/
    for (i = 0; i < m + 1; i++) { free(p_bar_numerators[i]); }
    free(p_bar_numerators);
    for (i = 0; i < m + 1; i++) { free(p_bar_denominators[i]); }
    free(p_bar_denominators);
    /*Build the tree*/
    T = build_tree_nj(distance_matrix, m);
    print_tree(T);
    /*Iterate through candidate loci to call variants*/
    candidate = malloc(sizeof(Candidate));
    candidate->m = m;
    rewind(f_candidates);
    read_candidate(f_candidates, candidate);
    /* TODO */
    delete_tree(T);
    /*Freeing memory, closing files*/
    for (i = 0; i < m + 1; i++) { free(distance_matrix[i]); }
    free(distance_matrix);
    free(p0); free(gp);
    free(P_sigma);
    free_log_factorials();
    free_log_binom();
    fclose(f_candidates);
    //FIXME uncomment after debugging
    //remove(gp->tmp_fname);

    return 0;
}

/* Stores candidate loci posterior distributions in a binary file */
int write_candidate(FILE* c_file, Locus* locus, long double P_0, long double* probs, int m) {
    int i;
    char is_valid;
    /*Length of sequence name */
    int seq_length = strlen(locus->sequence);
    fwrite(&seq_length, sizeof(int), 1, c_file);
    fwrite(locus->sequence, 1, seq_length, c_file);
    fputc('\0', c_file);
    fwrite(&(locus->position), sizeof(unsigned long), 1, c_file);
    fwrite(&P_0, sizeof(long double), 1, c_file);
    for (i = 0; i < m; i++) {
        is_valid = (locus->cells[i].read_count > 0) ? 1 : 0;
        fputc(is_valid, c_file);
        fwrite(probs + 3*i, sizeof(long double), 3, c_file);
    }
    return 0;
}

int read_candidate(FILE* c_file, Candidate* candidate) {
    int seq_length, i;
    fread(&seq_length, sizeof(int), 1, c_file);
    candidate->seq_name = malloc(seq_length + 1);
    fread(candidate->seq_name, 1, seq_length + 1, c_file);
    fread(&(candidate->pos), sizeof(unsigned long), 1, c_file);
    fread(&(candidate->P_sig0), sizeof(long double), 1, c_file);
    for (i = 0; i < m; i++) {
        valid[i] = fgetc(c_file);
        fread(probs + 3*i, sizeof(long double), 3, c_file);
        /*FIXME to candidate*/
    }
    
    free(sequence);
    return 0;
}

/*The expected number of alleles that differ between the two cells at a locus*/
long double pair_exp_difference(long double* posteriors, int a, int b) {
    long double to_sum[6];
    /*Genotypes differ by 2*/
    to_sum[0] = logl(2) + posteriors[3*a + 0] + posteriors[3*b + 2];
    to_sum[1] = logl(2) + posteriors[3*a + 2] + posteriors[3*b + 0];
    /*Genotypes differ by 1 */
    to_sum[2] = posteriors[3*a + 0] + posteriors[3*b + 1];
    to_sum[3] = posteriors[3*a + 1] + posteriors[3*b + 0];
    to_sum[4] = posteriors[3*a + 1] + posteriors[3*b + 2];
    to_sum[5] = posteriors[3*a + 2] + posteriors[3*b + 1];
    /*Numerator is the sum of these*/
    return LSE(to_sum, 6);
}


/*Returns number of candidates added from this batch*/
int update_candidates(Locus* loci_batch,
        int batch_size,
        prior_params_t* p,
        long double* sig_priors,
        FILE* candidates_file,
        long double** numr,
        int** denom) {
    int i, j, k, a, b, n_valid_cells;
    long double pair_exp_diff;
    int candidates_found = 0;
    long double* genotype_posteriors;
    int* valid_cells = malloc(p->m * sizeof(int));
    long double** cell_ls  = malloc(batch_size * sizeof(long double*));
    /*Posterior sigma distribution for each locus in batch*/
    long double** l_P_sig__D = malloc(batch_size * sizeof(long double*));
    long double*  locus_ls = malloc((2 * p->m + 1) * sizeof(long double));
    /*Iterate through loci for cell likelihoods and posterior sigma distributions*/
    for (i = 0; i < batch_size; i++) {
        cell_ls[i] = malloc(3 * p->m * sizeof(long double));
        for (j = 0; j < p->m; j++) {
            /*NB this function cleans reads of unknown and low quality, changing the cell struct*/
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
    /*Identify candidate loci write their genotype posteriors to tmp file
     * Also calculate P_bar numerators and denominators*/
    genotype_posteriors = malloc(3 * p->m * sizeof(long double));
    for (i = 0; i < batch_size; i++) {
        cell_posteriors(genotype_posteriors, l_P_sig__D[i], cell_ls[i], p->m);
        if (l_P_sig__D[i][0] <= p->thresh) {
            candidates_found++;
            write_candidate(candidates_file, &(loci_batch[i]), l_P_sig__D[i][0], genotype_posteriors, p->m);
        }
        /*P_bar calculations for building cell distances*/
        n_valid_cells = 0;
        /*See what cells are valid*/
        for (j = 0; j < p->m; j++) {
            if(loci_batch[i].cells[j].read_count > 0) {
                valid_cells[n_valid_cells++] = j;
            }
        }
        /*Iterate through pairs of valid cells*/
        for (j = 0; j < n_valid_cells; j++) {
            a = valid_cells[j];
            for (k = j+1; k < n_valid_cells; k++) {
                b = valid_cells[k];
                denom[a][b]++; denom[b][a]++;
                pair_exp_diff = pair_exp_difference(genotype_posteriors, a, b);
                numr[a][b] = LSE2(numr[b][a], pair_exp_diff);
                numr[b][a] = numr[a][b];
            }
            /*For each valid cell also compute distance with reference (outgroup)*/
            pair_exp_diff = LSE2(logl(2) + genotype_posteriors[3*a + 2], genotype_posteriors[3*a + 1]);
            numr[a][p->m] = numr[p->m][a] = LSE2(numr[a][p->m], pair_exp_diff);
            denom[a][p->m]++; denom[p->m][a]++; 
        }
        free(l_P_sig__D[i]);
        free(cell_ls[i]);
    }
    free(l_P_sig__D);
    free(cell_ls);
    free(valid_cells);
    free(genotype_posteriors);
    return candidates_found;
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
    double threshold = 0.5;
    strcpy(gp->tmp_fname, "/tmp/SCarborSNV_cand_tmp");
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
    /*TODO add options for threshold, tempfilename, amplification err and P_ADO*/
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
    p0->thresh      = log(threshold);
    p0->m           = m_cells;
}

