/**
 * SCarborSNV: Efficient Phylogeny-aware Single Nucleotide Variant Detection for Single Cells
 *
 * Copyright (C) 2019 Christopher Oldham
 *
 * This file is part of SCarborSNV.
 *
 * SCarborSNV is free software: you can redistribute it and/or modify
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SCarborSNV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SCarborSNV.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include "sequence_utils.h"

nuc_t decode_nucleotide(char encoded, nuc_t ref) {
    nuc_t ret = INVALID_NUC;
    switch(encoded) {
        case 'A' :
        case 'a' : ret = A;
        break;
        case 'C' :
        case 'c' : ret = C;
        break;
        case 'G' :
        case 'g' : ret = G;
        break;
        case 'T' :
        case 't' : ret = T;
        break;
        case '.' :
        case ',' : ret = ref;
        break;
        case '*' : ret = DELETED_NUC;
        break;
        case 'N' :
        case 'n' : ret = UNKNOWN_NUC;
        break;

    }
    return ret;
}

nuc_t decode_ref(char encoded) {
    nuc_t ret = INVALID_NUC;
    switch(encoded) {
        case 'a' :
        case 'A' : ret = A;
        break;
        case 'c' :
        case 'C' : ret = C;
        break;
        case 'g' :
        case 'G' : ret = G;
        break;
        case 't' :
        case 'T' : ret = T;
        break;
        case 'n' :
        case 'N' : ret = UNKNOWN_NUC;
        break;
    }
    return ret;
}

long double char2l_err(char c) {
    long double phred = (long double)(c - 33);
    return phred * phred_multiplier;
}

int clean_fill(int read_depth, nuc_t ref, char* raw_reads, char* raw_quals, nuc_t* processed_reads, long double* processed_quals) {
    int i, j, indel_length;
    char decimal_string[32];
    char* read_ptr = raw_reads;
    char* qual_ptr = raw_quals;
    if (read_depth == 0) return 0;
    for (i = 0; i < read_depth; i++) {
        /*Start segment + mapping quality */
        if (*read_ptr == '^') {
            read_ptr += 2;
            i--;
            continue;
        }
        /*End segment */
        else if (*read_ptr == '$') {
            read_ptr++;
            i--;
            continue;
        }
        /*Indel */
        else if (*read_ptr == '+' || *read_ptr == '-') {
            /*Number after +/- is indel length */
            read_ptr++;
            j = 0;
            while ('0' <= *read_ptr && *read_ptr <= '9') {
                decimal_string[j] = *read_ptr;
                read_ptr++;
                j++;
            }
            decimal_string[j] = '\0';
            indel_length = 0;
            indel_length = atoi(decimal_string);
            /*Indel reads are next, skipping to next non-indel */
            read_ptr += indel_length;
            /*TODO test this */
            i--;
            continue;
        }
        /*Normal read */
        else {
            processed_reads[i] = decode_nucleotide(*(read_ptr++), ref);
            processed_quals[i] = char2l_err(*(qual_ptr++));
        }
    } /*End for */
    if (*qual_ptr != '\0' || *read_ptr != '\0') {
        /*Could have characters beyond final valid read if indel or sequence markers */
        if (*read_ptr == '$' || *read_ptr == '^' || *read_ptr == '+' || *read_ptr == '-') return 0;
            return 1;
    }
    return 0;
}

nuc_t get_alt_allele(Locus* locus, nuc_t ref, int m) {
    int i, j, maxc;
    int counts[4] = {0, 0, 0, 0};
    for (i = 0; i < m; i++) {
        for (j = 0; j < locus->cells[i].read_count; j++) {
            switch (locus->cells[i].reads[j]) {
                case 0 : counts[0]++;
                         break;
                case 1 : counts[1]++;
                         break;
                case 2 : counts[2]++;
                         break;
                case 3 : counts[3]++;
                         break;
            }
        }
    }
    j = maxc = -1;
    for (i = 0; i < 4; i++) {
        if ((counts[i] > maxc) && ((int)ref != i)) {
            j = i;
            maxc = counts[i];
        }
    }

    return (nuc_t)j;
}

    

char base2char(nuc_t base) {
    switch(base)  {
        case 0 : return 'A';
        case 1 : return 'C';
        case 2 : return 'G';
        case 3 : return 'T';
        case 4 : return 'X';
        case 6 : return 'N';
    }
    return 'Z';
}
