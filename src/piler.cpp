#include <iostream>
#include <fstream>
#include <sstream>
#include "scarborsnv.h"
#include "piler.h"
#include "sequence_utils.h"

using namespace piler_module;
//TODO is order important? Or do we keep track of chrom and pos?

/*********************************
 * Piler methods
 * ******************************/

Piler::Piler(std::istream* instream, int n_threads, int batch_size, bool is_bams) {
    this->batch_size = batch_size;
    this->pileup_stream = instream;
    this->last_batch_id = -1;
    //TODO real queue size
    this->batch_queue.set_max_size(n_threads + 1);
    //TODO:
    std::string test;
    getline(*instream, test);
    std::cout << test << std::endl;
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
        Locus_reads* locus = new Locus_reads(chrom, pos, ref, n_cells);

        for (int i=0; i<n_cells; i++) {
            //TODO clean reads, convert phreds etc.
        }
    }
    return batch;

}

void Piler::fill_queue() {
    std::string line;
    //FIXME only one pileup stream, how can multiple threads work?
    //Maybe only multithread if bams. Only one thread can read mpileup stream at a time and thus
    //only one thread can make batches. Or we just keep track of locus positions, then it doesn't matter if they are in order.
    //FIXME when e.g. get_next_batch calls this, should skip unless it is the first worker doing it.
    while (this->batch_queue.size() < this->batch_queue.get_max_size()) {
        //Create a batch and add to queue


    }
    return;
}


Piler::~Piler() {
    if (this->pileup_stream != &std::cin) {
        std::ifstream* i =  static_cast<std::ifstream*>(this->pileup_stream);
        i->close();
    }
}


