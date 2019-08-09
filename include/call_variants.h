#ifndef CALL_VARIANTS_H_
#define CALL_VARIANTS_H_

#include "posteriors.h"

#define CALL_X (-1)
#define CALL_0 (0)
#define CALL_1 (1)
#define CALL_2 (2)

/*Writes a VCF line based on phylogeny aware posteriors and optionally initial posteriors.
 * To only include phylogeny aware calls, set zero_thresh to logl(0)
 * Phylogeny aware calls are made by maximum posterior genotype probability */
int call_to_VCF(FILE* vcf_file, Candidate* c, long double zero_thresh);

#endif
