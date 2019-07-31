#include <iostream>
#include <fstream>
#include "scarborsnv.h"
#include "piler.h"
#include "sequence_utils.h"

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
    this->cells = new Cell[n_cells];
    return;
}

void Locus::add_cell(int depth, std::string read_string, std::string qual_string) {
    //Cells should already be dynamically allocated by constructor
    //Read array in Cell allocated here?
    int cell_id = this->last_cell = this->last_cell +1;
    Cell* c = this->get_cell(cell_id);
    c->depth = depth;
    c->reads = new read[depth];
    sequence_utils::clean_fill(c->reads, depth, this->ref, read_string, qual_string);
    return;
}

Cell* Locus::get_cell(int i) {
    return (this->cells) + i;
}

std::string Locus::get_chrom() {
    return this->chrom;
}

unsigned int Locus::get_pos() {
    return this->pos;
}

nuc_t Locus::get_ref(){
    return this->ref;
}

Locus::~Locus() {
    for (int i=0; i<this->n_cells; i++) {
        delete[] this->cells[i].reads;
    }
    delete[] this->cells;
    return;
}

