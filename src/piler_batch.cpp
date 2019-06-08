#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
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
    return this->batch_id;

}

Locus* Batch::get_locus(int i) {
    return this->data[i];
}

void Batch::set_locus(int i, Locus* locus) {
    this->data[i] = locus;
}

void Batch::resize(unsigned int new_size) {
    Locus** new_data = new Locus*[new_size];
    //The Loci beyond have not yet been allocated, no need to delete just set new array to point at existing loci
    for (unsigned int i=0; i< new_size; i++) {
        new_data[i] = this->data[i];
    }
    delete[] this->data;
    this->data = new_data;
    this->batch_size = new_size;
    return;
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
    this->pileup_complete = false;
}

void Batch_Q::set_max_size(unsigned int max_size) {
    this->max_size = max_size;
    return;
}

Batch* Batch_Q::pop() {
    //TODO TESTME
    std::unique_lock<std::mutex> lock1(this->queue_mutex, std::defer_lock);
    Batch* batch;
    while (true) {
        lock1.lock();
        //Queue is empty but more on the way, wait
        if (this->batch_queue.empty() && !this->pileup_complete) {
            lock1.unlock();
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
        else if (this->batch_queue.empty() && this->pileup_complete) {
            lock1.unlock();
            batch = NULL;
            break;
        }
        //Queue not empty
        else {
            batch = this->batch_queue.front();
            this->batch_queue.pop();
            lock1.unlock();
            break;
        }
    }
    return batch;
}

void Batch_Q::push(Batch* b){
    std::lock_guard<std::mutex> lock1(this->queue_mutex);
    this->batch_queue.push(b);
}

unsigned int Batch_Q::size() {
    std::lock_guard<std::mutex> lock1(this->queue_mutex);
    return this->batch_queue.size();
}

unsigned int Batch_Q::get_max_size() {
    return this->max_size;
}

Batch_Q::~Batch_Q() {
    std::lock_guard<std::mutex> lock1(this->queue_mutex);
    Batch* batch;
    while (!this->batch_queue.empty()) {
        batch = this->batch_queue.front();
        this->batch_queue.pop();
        //FIXME this should rarely be used. Batches should already be deleted
        //delete batch;
    }
    return;
}
