#include <stdlib.h>
#include <string.h>
#include "pileup_reader.h"
#include "sequence_utils.h"

/*Maximum length of sequence name*/
#define L_SEQNAME (255)

char* dynamic_read(FILE* instream, char* read_buffer, size_t expected_size) {
    char c;
    size_t size = expected_size;
    char* buffer = read_buffer;
    int found_start = 0;
    int chars_read = 0;
    while (1) {
        c = fgetc(instream);
        if(!found_start && (c == '\t' || c == ' ')) {
            continue;
        }
        found_start = 1;
        if (c == '\t' || c == ' ') {
            buffer[chars_read] = '\0';
            return buffer;
        }
        buffer[chars_read] = c;
        chars_read += 1;
        if (chars_read + 4 >= size) { 
            buffer = realloc(buffer, size * 2);
            size *= 2;
        }
    }
}


int read_locus(FILE* instream, int m, Locus* locus) {
    int i;
    char* read_buffer;
    char* qual_buffer;
    char raw_ref;

    /*Name of sequence, first mpileup field*/
    locus->sequence = (char*)malloc(L_SEQNAME * sizeof(char));
    /*Other locus-wide data */
    if (fscanf(instream, "%s %lu %c" , locus->sequence, &(locus->position), &raw_ref) != 3) {
        free(locus->sequence);
        locus->ref_base = (feof(instream)) ? PILEEOF_NUC : INVALID_NUC;
        return 1;
    }
    locus->ref_base = decode_ref(raw_ref);

    /*Array of m Cell_locus structs */
    Cell_locus* cells = (Cell_locus*)malloc(m * sizeof(Cell_locus));
    for (i = 0; i < m; i++) {
        fscanf(instream, "%d", &(cells[i].read_count));
        cells[i].cell_position = i;
        /*Indels can cause read string to be much longer than depth*/
        read_buffer = (char*)malloc((cells[i].read_count + 32) * sizeof(char));
        read_buffer = dynamic_read(instream, read_buffer, cells[i].read_count + 32);
        /*Qual string should always only be the length of the read depth*/
        qual_buffer = (char*)malloc((cells[i].read_count + 32) * sizeof(char));
        fscanf(instream, "%s", qual_buffer);
        cells[i].reads = (nuc_t*)malloc((cells[i].read_count + 2) * sizeof(nuc_t));
        cells[i].quals = (long double*)calloc((cells[i].read_count + 2) , sizeof(long double));
        clean_fill(cells[i].read_count, locus->ref_base, read_buffer, qual_buffer, cells[i].reads, cells[i].quals);
        free(read_buffer);
        free(qual_buffer);
    }
    locus->cells = cells;

    return 0;

 }

int read_batch_loci(FILE* instream, Locus* loci, int n, int m) {
    int i;
    int status = 0;
    for (i = 0; i < n; i++) {
        status = read_locus(instream, m, &loci[i]);
        if (status != 0) {
            break;
        }
    }

    return i;
}

int delete_locus_contents(Locus* loci, int n, int m) {
    /*XXX only delete those have been filled.
     * Pass n as this number, received from reader*/
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            free(loci[i].cells[j].reads);
            free(loci[i].cells[j].quals);
        }
        free(loci[i].cells);
        free(loci[i].sequence);
    }
    return 0;
}

