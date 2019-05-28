#ifndef PILER_BATCH_H_
#define PILER_BATCH_H_

#include <vector>
#include <queue>
#include <mutex>
#include "scarborsnv.h"

namespace piler_module {

/* All the read and qual data from a single cell at a locus*/
class Locus_reads {
    public:
        Locus_reads(std::string chrom, unsigned int pos, nuc_t ref, int n_cells);
        void add_cell(int n_reads, read* cell_reads);
        ~Locus_reads();


    private:
        std::string chrom;
        unsigned int pos;
        int n_cells, last_cell;
        nuc_t ref;
        //read defined in scarborsnv.h
        read** reads;
};

/* A batch of 'batch_size' loci with the reads of each cell at those loci.
 * Note Batches get made and analyzed by single threads so don't need mutexes
 * */
class Batch {
    public:
        /* Allocates memory for the batch*/
        Batch(unsigned int batch_size, int id);
        /*Returns size of this batch even if smaller than standard batch size*/
        unsigned int get_batch_size();
        int get_batch_id();
        Locus_reads* get_locus(int i);
        void set_locus(int i, Locus_reads* locus);
        ~Batch();

    private:
        unsigned int batch_size;
        int batch_id;
        Locus_reads** data;
};

/*Thread safe queue of batches ready for pickup by workers */
class Batch_Q {
    public:
        Batch_Q();
        void set_max_size(unsigned int max_size);
        unsigned int get_max_size();
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
        //TODO initialize, use lock?
        std::mutex m_queue;
};


} //end namespace
#endif
