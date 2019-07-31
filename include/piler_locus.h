#ifndef PILER_LOCUS_H_
#define PILER_LOCUS_H_

#include <string>
#include "scarborsnv.h"

namespace piler_module {

/* A single cell at a locus */
struct Cell {
    int depth;
    //read defined in scarborsnv.h
    read* reads;
};

class Locus {
    public:
        Locus(std::string chrom, unsigned int pos, nuc_t ref, int n_cells);
        void add_cell(int depth, std::string read_string, std::string qual_string);
        Cell* get_cell(int i);
        std::string get_chrom();
        unsigned int get_pos();
        nuc_t get_ref();
        ~Locus();


    private:
        std::string chrom;
        unsigned int pos;
        nuc_t ref;
        int n_cells, last_cell;
        Cell* cells;
};

} //end namespace
#endif
