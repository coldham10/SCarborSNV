#include <string>
#include <stdexcept>
#include "sequence_utils.h"

//TODO clean string of inserts, etc

nuc_t sequence_utils::decode_nucleotide(char encoded, nuc_t ref) {
    char ret;
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
    }
    if (ret == INVALID_NUC) {
        throw std::invalid_argument("encoded reference (ACGT)");
    }
    return ret;
}
