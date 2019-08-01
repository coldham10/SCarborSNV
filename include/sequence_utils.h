#ifndef SEQUENCE_UTILS_H_
#define SEQUENCE_UTILS_H_

#include <tgmath.h>
#include "scarborsnv.h"


//Multiply a phred by this to get the log probability of error
static const long double phred_multiplier = -1/(10 * log10l(expl(1)));

nuc_t decode_ref(char encoded);

/* Converts raw pileup read and qual strings into nuc_t's and log qualities*/
int clean_fill(int read_depth, nuc_t ref, char* raw_reads, char* raw_quals, nuc_t* processed_reads, long double* processed_quals);

/* Converts a nuc_t base to a human readable character *
char base2char(nuc_t base);
*/

#endif
