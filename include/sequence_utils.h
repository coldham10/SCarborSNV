#ifndef SEQUENCE_UTILS_H_
#define SEQUENCE_UTILS_H_

#include <tgmath.h>
#include <string>
#include "scarborsnv.h"

namespace sequence_utils {

//Multiply a phred by this to get the log probability of error
static const long double phred_multiplier = -1/(10 * log10l(expl(1)));

/* Turns chars 'A' etc into the corresponding nuc_t*/
nuc_t decode_nucleotide(char encoded, nuc_t ref);
nuc_t decode_ref(char encoded);

/* Converts a pileup quality character to a log error */
long double char2l_err(char& c);

/* Converts raw pileup read and qual strings into an array of read structs */
void clean_fill(read* to_fill, int read_depth, nuc_t ref, std::string read_string, std::string qual_string);

/* Converts a nuc_t base to a human readable character */
char base2char(nuc_t base);

} //end namespace sequence_utils

#endif
