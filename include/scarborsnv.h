#ifndef SCARBORSNV_H_
#define SCARBORSNV_H_

typedef struct {
    unsigned int n_threads;
    unsigned int mp_isfile; //piling up is to be handled by SCarborSNV //TODO if many flags can do bit field
//    char mp_fname[256];
    std::string mp_fname;
} global_params;


#endif  //SCARBORSNV_H_
