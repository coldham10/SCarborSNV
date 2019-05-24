#ifndef SCARBORSNV_H_
#define SCARBORSNV_H_

static const char A = 0;
static const char C = 1;

typedef struct {
    unsigned int n_threads;
    bool mp_isfile; //TODO: piling up is to be handled by SCarborSNV
    std::string mp_fname;
} global_params_t;


#endif  //SCARBORSNV_H_
