#ifndef PILER_LOCUS_H_
#define PILER_LOCUS_H_

#include "scarborsnv.h"

namespace piler_module {

class Cell {

};

/* All the read and qual data from a single cell at a locus*/
class Locus {
    public:
        Locus(std::string chrom, unsigned int pos, nuc_t ref, int n_cells);
        void add_cell(int n_reads, read* cell_reads);
        read* get_cell(int i);
        void set_cell(int i, read* cell);
        ~Locus();


    private:
        std::string chrom;
        unsigned int pos;
        int n_cells, last_cell;
        nuc_t ref;
        //read defined in scarborsnv.h
        read** reads;
};

} //end namespace
#endif
