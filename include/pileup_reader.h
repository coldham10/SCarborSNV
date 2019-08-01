#ifndef PILEUP_READER_H_
#define PILEUP_READER_H_

#include "scarborsnv.h"

typedef struct {
    nuc_t* reads;
    long double* quals;
    int read_count;
    int cell_position;
} Cell_locus;

typedef struct {
    char* sequence;
    Cell_locus* cells;
    unsigned long position;
    nuc_t ref_base;
} Locus;

/*Read n loci from instream into the array starting at loci.
 * returns the number of loci successfully allocated and read*/
int read_batch_loci(FILE* instream, Locus* loci, int n, int m);

/*Deletes anything alloc'd by read_batch_loci. WARNING: still need to 
 * free the array of Locus separately*/
int delete_locus_contents(Locus* loci, int n, int m);

#endif /*define PILEUP_READER_H_ */
