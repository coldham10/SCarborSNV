#ifndef INFERENCE_H_
#define INFERENCE_H_

#include "tree.h"

/*Runs the upwards-downwards algorithm to prepare the tree in O(m) then computes genotype probabilities by DP*/
int infer_from_phylogeny(Node* T, long double* posteriors, long double* result, long double P_SNV, long double P_H);



#endif
