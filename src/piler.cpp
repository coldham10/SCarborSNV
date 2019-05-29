#include <iostream>
#include <fstream>
#include <sstream>
#include "scarborsnv.h"
#include "piler.h"
#include "piler_locus.h"
#include "sequence_utils.h"

using namespace piler_module;
//TODO is order important? Or do we keep track of chrom and pos?

/*********************************
 * Public Piler methods
 * ******************************/

Piler::Piler(std::istream* instream, int n_threads, int batch_size, bool is_bams) {
    this->batch_size = batch_size;
    //FIXME:
    this->batch_size = 2;
    this->pileup_stream = instream;
    this->last_batch_id = -1;
    this->batch_queue.set_max_size(n_threads + 1);
    Batch* test = this->make_batch();
    Locus* t2 = test->get_locus(0);
    Cell* t3 = t2->get_cell(1);
    printf("base is %d\n", static_cast<int>(t3->reads[1].base));
    //TODO
    delete test;
}

unsigned int Piler::get_batch_size() {
    return this->batch_size;
}

Batch* Piler::get_next_batch() {
    if (this->batch_queue.size() < this->min_queue_size) {
        //TODO: split off a worker to do this!:
        this->fill_queue();
    }
    return this->batch_queue.pop();
}

Piler::~Piler() {
    if (this->pileup_stream != &std::cin) {
        std::ifstream* i =  static_cast<std::ifstream*>(this->pileup_stream);
        i->close();
    }
}

/*********************************
 * Private Piler methods
 * ******************************/

Batch* Piler::make_batch() {
    Batch* batch = new Batch(this->batch_size, ++(this->last_batch_id));
    //Iterate over loci/pileup lines
    for (unsigned int i=0; i < this->batch_size; i++) {
        std::string line;
        getline(*(this->pileup_stream), line);

        std::vector<std::string> tokens;
        std::stringstream linestream(line);
        std::string token;
        //Split mpileup line into tokens
        while (getline(linestream, token, '\t')) {
            tokens.push_back(token);
        }

        std::string chrom = tokens[0];
        //TODO test on overflowing ints
        unsigned int pos = static_cast<unsigned int>(std::stoll(tokens[1]));
        //Three fields per cell + 3 fields for whole locus
        int n_cells = (tokens.size()/3) -1;
        nuc_t ref = sequence_utils::decode_nucleotide(tokens[2].front());
        Locus* locus = new Locus(chrom, pos, ref, n_cells);

        int depth;
        std::string read_string, qual_string;
        for (int j=0; j<n_cells; j++) {
            depth = stoi(tokens[3 + 3*j]);
            read_string = tokens[3 + 3*j + 1];
            qual_string = tokens[3 + 3*j + 2];

            locus->add_cell(depth, read_string, qual_string);
            //sequence_utils::clean_fill(locus->get_cell(j), depth, ref, read_string, qual_string);
        }
        batch->set_locus(i, locus);
    }
    return batch;

}

void Piler::fill_queue() {
    std::string line;
    while (this->batch_queue.size() < this->batch_queue.get_max_size()) {
        //Create a batch and add to queue


    }
    return;
}




