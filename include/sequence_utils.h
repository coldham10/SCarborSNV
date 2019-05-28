#ifndef SEQUENCE_UTILS_H_
#define SEQUENCE_UTILS_H_

#include "scarborsnv.h"

namespace sequence_utils {

//TODO clean string of inserts, etc

/* Turns chars 'A' etc into the corresponding nuc_t*/
nuc_t decode_nucleotide(char encoded, nuc_t ref=INVALID_NUC);

} //end namespace sequence_utils

#endif
