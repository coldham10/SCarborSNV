#ifndef PILER_H_
#define PILER_H_

#include <iostream>

class Piler {
    public:
        Piler(std::istream& instream, bool is_bams=false);
};

/*Eventually three possible options:
  1. Read from pileup file
  2. Read pileup from piped in samtools mpileup
  3. Use htslib/samtools c api to take list of bams & pile 
*/


#endif
