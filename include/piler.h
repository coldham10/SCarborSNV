#ifndef PILER_H_
#define PILER_H_

#include <iostream>
#include <vector>

struct Cell_reads {
    std::vector<read> reads;
    int n_reads;
    nuc_t ref;
};

class Batch {
    public:
        Batch(int batch_size);
        ~Batch();
        int get_batch_size();
        Cell_reads get_cell(int i);

    private:
        int batch_size;
        Cell_reads* data;
};


class Piler {
    public:
        /* Sets the read stream and starts working to populate batch list
         * */ 
        Piler(std::istream* instream, int batch_size=1000, bool is_bams=false);
        ~Piler();
        //TODO destructor that closes file if not cin

        int get_m_cells();
        Batch get_next_batch();

    private:
        std::istream* pileup_stream;
        int m_cells, batch_size;
};

/*Eventually three possible options:
  1. Read from pileup file
  2. Read pileup from piped in samtools mpileup
  3. Use htslib/samtools c api to take list of bams & pile 
*/

//IDEA: have a list of batches ready (n=O(2*number of threads)). When initialized or batch requested, spawn a worker that works to repopulate
//the batches until there are enough, then dies. When a main worker requests a new batch, should have a whole batch ready to go already piled.
//?If not, dedicate more workers to piling? Will have to work on this...

#endif
