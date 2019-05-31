#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include "scarborsnv.h"
#include "piler.h"
#include "piler_locus.h"
#include "sequence_utils.h"

using namespace piler_module;
/*********************************
 * Public Piler methods
 * ******************************/

//FIXME does not work for files (spurious blank lines)
Piler::Piler(std::istream* instream, bool is_stdin, int n_cells, int n_threads, int batch_size, bool is_bams) {
    this->batch_size = batch_size;
    //XXX
    this->batch_size = 10;
    this->pileup_stream = instream;
    this->is_stdin = is_stdin;
    this->n_cells = n_cells;
    this->last_batch_id = -1;
    this->pileup_complete = false; //this->batch_queue.pileup_complete;
    this->batch_queue.set_max_size(n_threads + 1);
    
    this->try_fill_queue();
}

unsigned int Piler::get_batch_size() {
    return this->batch_size;
}

int Piler::get_n_batches() {
    return this->batch_queue.size();
}

Batch* Piler::get_next_batch() {
    if (this->batch_queue.size() < this->batch_queue.get_max_size()) {
        this->try_fill_queue();
    }
    return (this->batch_queue.size() == 0) ? NULL : this->batch_queue.pop();
}

Piler::~Piler() {
    if (!this->is_stdin) {
        std::ifstream* i =  static_cast<std::ifstream*>(this->pileup_stream);
        i->close();
    }
}

/*********************************
 * Private Piler methods
 * ******************************/

bool Piler::can_fill() {
    std::lock_guard<std::mutex> lock(this->f_active_mutex);
    if (this->f_active || this->pileup_complete) {
        return false;
    }
    else {
        this->f_active = true;
        //Caller can begin filling queue
        return true;
    }
}

void Piler::filler_done() {
    std::lock_guard<std::mutex> lock(this->f_active_mutex);
    this->f_active = false;
}


void Piler::try_fill_queue() {
    //NB checking can_fill() can change f_active
    if (this->batch_queue.size() < this->batch_queue.get_max_size() && this->can_fill()) {
        //FIXME should be done with a separate thread. Detached or somehow join()
        //this->current_filler = new std::thread(&Piler::fill_queue, this);
        this->fill_queue();

    }
    return;
}

void Piler::fill_queue() {
    //Only callable from try_fill_queue
    while (this->batch_queue.size() < this->batch_queue.get_max_size() && !this->pileup_complete) {
        Batch* b = this->make_batch();
        this->batch_queue.push(b);
    }
    this->filler_done();
    return;
}



Batch* Piler::make_batch() {
    Batch* batch = new Batch(this->batch_size, ++(this->last_batch_id));
    //Iterate over loci/pileup lines
    for (unsigned int i=0; i < this->batch_size; i++) {
        std::string line;
        getline(*(this->pileup_stream), line);
        //std::cout << line << std::endl;

        if (line.empty()) {
            //FIXME If is bams and no piles ready yet, would getline return empty?
            //We do not want to resize the batch in this case
            //TODO indicate to Piler the EOF
            batch->resize(i+1);
            this->pileup_complete = true;
            //TODO send message to batch_Q
            return batch;
        }

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
        nuc_t ref = sequence_utils::decode_nucleotide(tokens[2].front());
        Locus* locus = new Locus(chrom, pos, ref, this->n_cells);

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


