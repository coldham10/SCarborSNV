#include <iostream>
#include <fstream>
#include "scarborsnv.h"
#include "piler.h"

using namespace piler_module;

/*********************************
 * Piler methods
 * ******************************/

Piler::Piler(std::istream* instream, int n_threads, int batch_size, bool is_bams) {
    this->batch_size = batch_size;
    this->pileup_stream = instream;
    this->last_batch_id = -1;
    //TODO real queue size
    this->batch_queue.set_max_size(5);
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
    //TODO
    return;
}


Piler::~Piler() {
    if (this->pileup_stream != &std::cin) {
        std::ifstream* i =  static_cast<std::ifstream*>(this->pileup_stream);
        i->close();
    }
}


