#include <iostream>
#include <fstream>
#include "scarborsnv.h"
#include "piler_batch.h"
#include "piler_locus.h"
#include "piler.h"

using namespace piler_module;
/*********************************
 * Batch methods
 * ******************************/

Batch::Batch(unsigned int batch_size, int batch_id) {
    this->batch_id = batch_id;
    this->batch_size = batch_size;
    this->data = new Locus*[batch_size]; 

    return;
}

unsigned int Batch::get_batch_size() {
    return this->batch_size;
}

int Batch::get_batch_id() {
    //TODO
    return 0;

}

Locus* Batch::get_locus(int i) {
    return this->data[i];
}

void Batch::set_locus(int i, Locus* locus) {
    this->data[i] = locus;
}

Batch::~Batch() {
    for (unsigned int i=0; i<this->batch_size; i++) {
        //Locus*
        delete this->data[i];
    }
    delete[] this->data;

    return;
}

/*********************************
 * Batch_Q methods
 * ******************************/

Batch_Q::Batch_Q() {
    this->max_size = -1;
    this->read_complete = false;
    this->n = 0;
    //std::queue<Batch*> data;
    //TODO initialize
    //std::mutex m_queue;
    //TODO !!!
}

void Batch_Q::set_max_size(unsigned int max_size) {
    this->max_size = max_size;
    return;
}

Batch* Batch_Q::pop() {
    //TODO !!!
        /*Returns oldest batch and deletes from queue
         * If empty returns null pointer TODO:(?)*/
        //TODO: mutex protected. If no batch available yet(and !read_complete), what to do? return null? Wait? reader only takes mutex when pushing complete Batch. 
    return (Batch*) NULL;
}

void Batch_Q::push(Batch* b){
    //TODO
}

unsigned int Batch_Q::size() {
    //TODO mutex lock;
    return this->n;
}

unsigned int Batch_Q::get_max_size() {
    return this->max_size;
}

Batch_Q::~Batch_Q() {
    //TODO
}
