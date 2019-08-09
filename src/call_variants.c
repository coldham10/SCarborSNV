#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "call_variants.h"
#include "sequence_utils.h"

#define MAX_CALL_QUAL (100)




int call_to_VCF(FILE* vcf_file, Candidate* c, long double zero_thresh) {
    int i, j, max_j, n_variants, n_valid;
    long double maxv, qual;
    int m = c->m;
    char call_str[16];
    char* calls = malloc(m);
    for (i = 0; i < m; i++) { calls[i] = CALL_X; }
    n_variants = 0;
    for (i = 0; i < m; i++) {
        max_j = -1;
        maxv = -INFINITY;
        for (j = 0; j < 3; j++) {
            if (c->phylo_posteriors[3*i + j] >= maxv) {
                maxv = c->phylo_posteriors[3*i + j];
                max_j = j;
            }
        }
        if (max_j > 0) {
            calls[i] = (char)max_j;
            n_variants++;
        }
        else if (c->simple_posteriors[3*i] < zero_thresh && c->valid_cells[i]) {
            calls[i] = CALL_1;
            n_variants++;
        }
    }
    if (n_variants == 0) { return 1; }
    n_valid = 0;
    for (i = 0; i < m; i++) {
        n_valid += c->valid_cells[i];
    }
    qual = roundl(c->P_0 / phred_multiplier);
    qual = (qual > MAX_CALL_QUAL) ? MAX_CALL_QUAL : qual;
    fprintf(vcf_file, "%s\t%lu\t.\t%c\t%c\t%d\tGT", c->seq_name, c->pos, base2char(c->ref), base2char(c->alt), (int)qual);
    /*TODO add other fields such as depth, possibly distinguish between phylo calls and pure probability catches*/
    for (i = 0; i < m; i++) {
        switch(calls[i]) {
            case CALL_X : strcpy(call_str, "./.");
                          break;
            case CALL_0 : strcpy(call_str, "0/0");
                          break;
            case CALL_1 : strcpy(call_str, "0/1");
                          break;
            case CALL_2 : strcpy(call_str, "1/1");
                          break;
        }
        fprintf(vcf_file, "\t%s", call_str);
    }
    free(calls);
    fprintf(vcf_file, "\n");
    return 0;
}
