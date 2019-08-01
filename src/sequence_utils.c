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
/*void sequence_utils::clean_fill(read* to_fill, int read_depth, nuc_t ref, std::string read_string, std::string qual_string) { */
    int i, j, indel_length;
    char decimal_string[32];
    char* read_ptr = raw_reads;
    char* qual_ptr = raw_quals;
    if (read_depth == 0) return 0;
    for (i=0; i<read_depth; i++) {
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
            
/*
char sequence_utils::base2char(nuc_t base) {
    switch(base)  {
        case A : return 'A';
        case C : return 'C';
        case G : return 'G';
        case T : return 'T';
        case UNKNOWN_NUC : return 'N';
        case INVALID_NUC : return 'X';
    }
    return 'Z';
}
*/
