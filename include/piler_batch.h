#ifndef PILER_BATCH_H_
#define PILER_BATCH_H_

#include <iostream>
#include <vector>
#include <queue>
#include <mutex>

namespace piler_module {

/* All the read and qual data from a single cell at a locus*/
//TODO when are these allocated and deallocated? Make Locus_reads a class?
struct Cell_reads {
    std::vector<read> reads;
    int n_reads;
    nuc_t ref;
};

typedef std::vector<Cell_reads> Locus_reads;

/* A batch of 'batch_size' loci with the reads of each cell at those loci.
 * Note Batches get made and analyzed by single threads so don't need mutexes
 * */
class Batch {
    public:
        /* Allocates memory for the batch*/
        Batch(int batch_size, int id);
        /*Returns size of this batch even if smaller than standard batch size*/
        unsigned int get_batch_size();
        int get_batch_id();
        Locus_reads* get_locus(int i);
        ~Batch();

    private:
        int batch_size, batch_id;
        Locus_reads** data;
};

/*Thread safe queue of batches ready for pickup by workers */
class Batch_Q {
    public:
        Batch_Q();
        void set_max_size(int max_size);
        /*Returns oldest batch and deletes from queue
         * If empty returns null pointer TODO:(?)*/
        //TODO: mutex protected. If no batch available yet(and !read_complete), what to do? return null? Wait? reader only takes mutex when pushing complete Batch. 
        Batch* pop();
        void push(Batch* b);
        unsigned int size();
        ~Batch_Q();

    private:
        unsigned int max_size, n;
        bool read_complete;
        std::queue<Batch*> batch_queue;
        //TODO initialize
        std::mutex m_queue;
};


} //end namespace
#endif
