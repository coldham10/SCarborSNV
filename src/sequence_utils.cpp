#include <stdexcept>
#include "sequence_utils.h"

nuc_t sequence_utils::decode_nucleotide(char encoded, nuc_t ref) {
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
    if (ret == INVALID_NUC) {
        throw std::invalid_argument("Encoded reference (ACGT). Either no reference was passed or the reference was invalid.");
    }
    return ret;
}

nuc_t sequence_utils::decode_ref(char encoded) {
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
    if (ret == INVALID_NUC) {
        throw std::invalid_argument("Encoded reference must be (ACGT). Got: "+ encoded); 
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
    for (int i=0; i<read_depth; i++) {
        //Start segment + mapping quality
        if (*read_it == '^') {
            read_it += 2;
            i--;
            continue;
        }
        //End segment
        else if (*read_it == '$') {
            read_it++;
            i--;
            continue;
        }
        //Indel
        else if (*read_it == '+' || *read_it == '-') {
            //Number after +/- is indel length
            std::string num_str = "";
            read_it++;
            while ('0' <= *read_it && *read_it <= '9') {
                num_str.append(1, *read_it);
                read_it++;
            }
            int n = 0;
            try {
                n = stoi(num_str);
            }
            catch (...) {
                throw std::runtime_error("Expected number after indel, got: " + num_str + " (" + std::to_string(num_str.length()) + ")\n");
            }
            //Indel reads are next, skipping to next non-indel
            read_it += n;
            //TODO test this
            i--;
            continue;

        }
        else {
            to_fill[i].base     = sequence_utils::decode_nucleotide(*(read_it++), ref);
            to_fill[i].l_err    = sequence_utils::char2l_err(*(qual_it++));
        }
    }
    if (qual_it != qual_string.end() || read_it != read_string.end()) {
        //Could have characters beyond final valid read if indel or sequence markers
        if (*read_it == '$' || *read_it == '^' || *read_it == '+' || *read_it == '-') return;
        std::string err = "Read string, qual string and depth don't all agree\n";
        err.append("Reads: " + read_string + "\n");
        err.append("Quals: " + qual_string + "\n");
        throw std::runtime_error(err);
    }
    return;
}
            
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

