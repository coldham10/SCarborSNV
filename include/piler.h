#ifndef PILER_H_
#define PILER_H_

#include <iostream>
#include <vector>
#include <queue>
#include <mutex>

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



class Piler {
    public:
        /* Sets the read stream and starts working to populate batch list
         * */ 
        Piler(std::istream* instream, int batch_size=1000, bool is_bams=false);
        unsigned int get_batch_size();
        Batch* get_next_batch();
        ~Piler();

    private:
        std::istream* pileup_stream;
        int m_cells;
        unsigned int min_queue_size, max_queue_size, last_batch_id, batch_size;
        Batch_Q batch_queue;
};

/*Eventually three possible options:
  1. Read from pileup file
  2. Read pileup from piped in samtools mpileup
  3. Use htslib/samtools c api to take list of bams & pile 
*/

//IDEA: have a list of batches ready (n=O(2*number of threads)). When initialized or batch requested, spawn a worker that works to repopulate
//the batches until there are enough, then dies. When a main worker requests a new batch, should have a whole batch ready to go already piled.
//?If not, dedicate more workers to piling? Will have to work on this...
//
//ALSO: in general, how to limit number of threads. If some thread (or many in case of in-house piling) are running in the piler, some in the other workers, how to keep total at t?
//Some sort of shared n_current threads variable, or even a thread manager class with mutex protected modifiers for this number and methods that can spawn new workers if some available. NOTE: this is for WAAAAY later on

#endif
