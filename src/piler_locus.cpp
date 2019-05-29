#include <iostream>
#include <fstream>
#include "scarborsnv.h"
#include "piler.h"

using namespace piler_module;
/*********************************
 * Locus methods
 * ******************************/
Locus::Locus(std::string chrom, unsigned int pos, nuc_t ref, int n_cells) {
    this->chrom = chrom;
    this->pos = pos;
    this->ref = ref;
    this->n_cells = n_cells;
    this->last_cell = -1;
    this->reads = new read*[n_cells];
    return;
}

void Locus::add_cell(int n_reads, read* cell_reads) {
    //cell_reads should already be dynamically allocated by here
    int pos = this->last_cell = this->last_cell +1;
    this->reads[pos] = cell_reads;
    return;
}

read* Locus::get_cell(int i) {
    return this->reads[i];
}

void Locus::set_cell(int i, read* cell) {
    this->reads[i] = cell;
    return;
}

Locus::~Locus() {
    for (int i=0; i<this->n_cells; i++) {
        delete[] this->reads[i];
    }
    delete[] this->reads;
    return;
}

