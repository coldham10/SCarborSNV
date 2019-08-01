#ifndef SCARBORSNV_H_
#define SCARBORSNV_H_

typedef char nuc_t;

static const nuc_t A = 0;
static const nuc_t C = 1;
static const nuc_t G = 2;
static const nuc_t T = 3;
static const nuc_t INVALID_NUC = 4;
static const nuc_t DELETED_NUC = 5;
static const nuc_t UNKNOWN_NUC = 6; //For 'N' bases in pileup data
static const nuc_t PILEEOF_NUC = 7;

typedef struct {
    unsigned int n_threads;
    int m; //Number of cells
    int mp_isfile; //TODO: piling up is to be handled by SCarborSNV
    char mp_fname[255];
} global_params_t;

#endif  //SCARBORSNV_H_
