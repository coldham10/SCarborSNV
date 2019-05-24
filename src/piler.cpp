#include <iostream>
#include "scarborsnv.h"
#include "piler.h"

Piler::Piler(std::istream& instream, bool is_bams) { //TODO: need ampersand after istream?
    std::string test;
    getline(instream, test);
    std::cout << test << std::endl;
}

//TODO read_data should be some sort of piler init, then have a get_batch function. Maybe create and return a piler class with a get batch method?

