#ifndef PILER_H_
#define PILER_H_

#include <iostream>
#include <mutex>
#include "piler_batch.h"

namespace piler_module {

class Piler {
    public:
        /* Sets the read stream and starts working to populate batch list
         * */ 
        Piler(std::istream* instream, bool is_stdin, int n_cells, int n_threads, int batch_size=1000, bool is_bams=false);
        unsigned int get_batch_size();
        int get_n_batches();
        /*Returns null if no batches left*/
        Batch* get_next_batch();
        ~Piler();

    private:
        //Private members
        std::istream* pileup_stream;
        int n_cells;
        unsigned int max_queue_size, last_batch_id, batch_size;
        Batch_Q batch_queue;
        bool is_stdin;
        //TODO IDEA: use << instead of getline, maybe solve crashing and fstream issues...
        bool f_active; //Is a filler making batches
        bool pileup_complete; //The pileup stream has ended
        std::mutex f_active_mutex;

        //Private methods
        /*Checks against f_active and pileup_done to see if the current thread can start filling batches.
         * If sucsessful (f_active & pileup_done were false) it sets f_active to true and returns true.*/
        bool can_fill();
        //called when filling complete due to full queue;
        void filler_done();
        /*Spawns a worker and begins filling queue, if possible*/
        void try_fill_queue();
        /*Continues to add batches to the queue until full */
        void fill_queue();
        Locus* peek_line();
        Batch* make_batch();
};

/*Eventually three possible options:
  1. Read from pileup file
  2. Read pileup from piped in samtools mpileup
  3. Use htslib/samtools c api to take list of bams & pile 
*/

} //end namespace
#endif
