#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <getopt.h>
#include "scarborsnv.h"
#include "piler.h"
#include "sigma_priors.h"
#include "sequence_utils.h"

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv);

int main(int argc, char** argv) {
    prior_params_t p0{};
    global_params_t gp{};
    init_params(&gp, &p0, argc, argv);

    std::vector<double> priors = log_sigma_priors(&p0);
    std::ofstream pout;
    pout.open("priors.n");
    pout << "m: " << p0.m << std::endl;
    for (int i=0; i<2*(p0.m)+1; i++) {
        pout << priors[i] << std::endl;
    }
    pout.close();

    /*
    piler_module::Piler* piler;
    if (gp.mp_isfile) {
        //Read from file
        std::ifstream ifs(gp.mp_fname, std::ifstream::in);
        piler = new piler_module::Piler(&ifs, false, gp.m, gp.n_threads, 10, false);
    }
    else {
        //Read from stdin
        piler = new piler_module::Piler(&std::cin, true, gp.m, gp.n_threads, 10, false);
    }
    piler_module::Batch* batch = piler->get_next_batch();

    //Writing another pileup to see if we got everything
    std::ofstream mirr;
    mirr.open("mirror.pileup");
    while ( batch != NULL) {
        //std::cout << "got batch, " << piler->get_n_batches() << " left." << std::endl;
        for (int i=0; i<batch->get_batch_size(); i++) {
            piler_module::Locus* loc = batch->get_locus(i);
            mirr << loc->get_chrom() << "\t" << loc->get_pos() << "\t" << sequence_utils::base2char(loc->get_ref());
            for (int j=0; j<gp.m; j++) {
                piler_module::Cell* c = loc->get_cell(j);
                mirr << "\t";
                for (int k=0; k<c->depth; k++) {
                    mirr << sequence_utils::base2char(c->reads[k].base);
                }
                mirr << "\t";
                for (int k=0; k<c->depth; k++) {
                    mirr << (char) static_cast<int>((c->reads[k].l_err / sequence_utils::phred_multiplier )+33);
                }
            }
        }
        mirr << "\n";
        delete batch;
        batch = piler->get_next_batch();
    }



    delete piler;
    */
    return 0;
}

void init_params(global_params_t* gp, prior_params_t* p0, int argc, char** argv) {
    //Default values
    unsigned int n_threads = 1;
    int m_cells = 1;
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
                //If pileup comes from file, -p pfile.pileup or --pileup-file pfile.xx
                mp_isfile = true;
                mp_fname = optarg;
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
           
    //Populate global parameter struct
    gp->n_threads   = n_threads;
    gp->mp_isfile   = mp_isfile;
    gp->mp_fname    = mp_fname;
    gp->m           = m_cells;

    //Convert prior params to log space, populate intial parameters struct
    p0->l_lambda    = std::log(lambda0);
    p0->l_mu        = std::log(mu0);
    p0->l_P_H       = std::log(P_H0);
    p0->l_P_clonal  = std::log(P_clonal0);
    p0->m           = m_cells;
}

