#ifndef SEQUENCE_UTILS_H_
#define SEQUENCE_UTILS_H_

#include <math.h>
#include "scarborsnv.h"
#include "pileup_reader.h"


/*Multiply a phred by this to get the log probability of error. */
/* ln(P) = -(Q/10)*ln(10) */
static const long double phred_multiplier = -0.2302585092994;

nuc_t decode_ref(char encoded);

/* Converts raw pileup read and qual strings into nuc_t's and log qualities */
int clean_fill(int read_depth, nuc_t ref, char* raw_reads, char* raw_quals, nuc_t* processed_reads, long double* processed_quals);

/*Based on maximum count, gets the alternate allele*/
nuc_t get_alt_allele(Locus* locus, nuc_t ref, int m);

/*Converts a nuc_t base to a human readable character */
char base2char(nuc_t base);


#endif
