#ifndef SCARBORSNV_H_
#define SCARBORSNV_H_

typedef char nuc_t;
typedef struct {
    long double l_err;
    nuc_t base;
} read;

static const nuc_t A = 0;
static const nuc_t C = 1;
static const nuc_t G = 2;
static const nuc_t T = 3;

typedef struct {
    unsigned int n_threads;
    bool mp_isfile; //TODO: piling up is to be handled by SCarborSNV
    std::string mp_fname;
} global_params_t;


#endif  //SCARBORSNV_H_
