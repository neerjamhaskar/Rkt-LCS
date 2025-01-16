//
// Created by hasibih on 1/16/25.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void generate_fasta_file(const char *filename, int num_reads, int read_length) {
    // Open the file for writing
    FILE *fasta_file = fopen(filename, "w");
    if (fasta_file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Seed the random number generator
    srand((unsigned int)time(NULL));

    // DNA alphabet
    const char dna_alphabet[] = "ATGC";

    for (int i = 1; i <= num_reads; i++) {
        // Write the header line
        fprintf(fasta_file, ">read_%d\n", i);

        // Generate a random DNA sequence
        for (int j = 0; j < read_length; j++) {
            char base = dna_alphabet[rand() % 4]; // Pick a random base from ATGC
            fputc(base, fasta_file);
        }

        // End the sequence line
        fputc('\n', fasta_file);
    }

    // Close the file
    fclose(fasta_file);
}

int main() {
    // Parameters for FASTA generation
    const char *filename = "fake_reads.fasta";
    int num_reads = 500;
    int read_length = 20;

    // Generate the FASTA file
    generate_fasta_file(filename, num_reads, read_length);

    printf("FASTA file '%s' generated with %d reads of length %d.\n", filename, num_reads, read_length);

    return 0;
}