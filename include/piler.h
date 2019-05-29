#ifndef PILER_H_
#define PILER_H_

#include <iostream>
#include "piler_batch.h"

namespace piler_module {

class Piler {
    public:
        /* Sets the read stream and starts working to populate batch list
         * */ 
        Piler(std::istream* instream, int n_threads, int batch_size=1000, bool is_bams=false);
        unsigned int get_batch_size();
        /*Returns null if no batches left*/
        Batch* get_next_batch();
        ~Piler();

    private:
        //Private members
        std::istream* pileup_stream;
        int m_cells;
        unsigned int min_queue_size, max_queue_size, last_batch_id, batch_size;
        Batch_Q batch_queue;
        bool f_active;

        //Private methods
        /*Checks against f_active to see if the current thread can start filling batches.
         * If sucsessful (f_active was false) it sets f_active to true and returns true.
         * Else returns false. TODO */
        bool can_fill();
        /*TODO spawns a worker and begins filling queue if possible*/
        void try_fill_queue();
        /* Continues to add batches to the queue until full */
        void fill_queue();
        Batch* make_batch();
        //TODO is this necessary? Just call delete?
        void delete_batch(Batch* batch);
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

} //end namespace
#endif
