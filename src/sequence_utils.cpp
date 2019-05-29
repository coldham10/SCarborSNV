#include <stdexcept>
#include "sequence_utils.h"

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
        throw std::invalid_argument("Encoded reference (ACGT). Either no reference was passed or the reference was invalid.");
    }
    return ret;
}

long double sequence_utils::char2l_err(char& c) {
    long double phred = static_cast<long double>(c - 33);
    return phred * sequence_utils::phred_multiplier;
}

void sequence_utils::clean_fill(read* to_fill, int read_depth, nuc_t ref, std::string read_string, std::string qual_string) {
    if (read_depth == 0) return;
    std::string::iterator read_it = read_string.begin();
    std::string::iterator qual_it = qual_string.begin();
    //int i = 0;
    //while (read_it != read_string.end() && qual_it != qual_string.end()) {
    for (int i=0; i<read_depth; i++) {
        //Start segment + mapping quality
        if (*read_it == '^') {
            read_it += 2;
            i--;
            continue;
        }
        //End segment
        else if (*read_it == '&') {
            read_it ++;
            i--;
            continue;
        }
        else {
            to_fill[i].base     = sequence_utils::decode_nucleotide(*(read_it++), ref);
            to_fill[i].l_err    = sequence_utils::char2l_err(*(qual_it++));
        }
    }
    if (read_it != read_string.end() || qual_it != qual_string.end()) {
        throw std::runtime_error("Read string, qual string and depth don't all agree");
    }
    return;
}
            


