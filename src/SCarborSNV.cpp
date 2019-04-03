#include <iostream>
#include "sam.h"

int main(void) {
    samFile* fp = sam_open("test.sam", "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    char* t = hdr->text;
    std::cout << t << std::endl;
    sam_close(fp);
    return 0;
}
