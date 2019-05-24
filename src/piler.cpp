#include <iostream>
#include <fstream>
#include "scarborsnv.h"
#include "piler.h"

Piler::Piler(std::istream* instream, int batch_size, bool is_bams) {
    //TODO: if_bams?
    this->batch_size = batch_size;
    this->pileup_stream = instream;
    std::string test;
    getline(*instream, test);
    std::cout << test << std::endl;
}

Piler::~Piler() {
    if (this->pileup_stream != &std::cin) {
        std::ifstream* i =  static_cast<std::ifstream*>(this->pileup_stream);
        i->close();
    }
    std::string test;
    getline(std::cin, test);
    std::cout << test << std::endl;
}

Batch::Batch(int batch_size) {
    Cell_reads* data = new Cell_reads[batch_size]; 
    for (int i=0; i<batch_size; i++) {
        data[i] = Cell_reads();
    }
}
