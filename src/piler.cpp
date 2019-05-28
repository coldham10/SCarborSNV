#include <iostream>
#include <fstream>
#include "scarborsnv.h"
#include "piler.h"

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
        //TODO: handle queue refillling here as necessary if queue below min size
    }
    return this->batch_queue.pop();
}

void Piler::fill_queue() {
    std::string line;
    //FIXME only one pileup stream, how can multiple threads work?
    //Maybe only multithread if bams. Only one thread can read mpileup stream at a time and thus
    //only one thread can make batches. Or we just keep track of locus positions, then it doesn't matter if they are in order.
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


