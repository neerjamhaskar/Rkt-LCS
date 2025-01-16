////  ./dependencies/bin/mpicc -o Flouri_OpenMPI.o Flouri_OpenMPI.c
///// ./dependencies/bin/mpirun -np 4 ./Flouri_OpenMPI.o test.FASTA 3
#include "dependencies/include/mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "flouri_lcp_table.h" // Changed to include table.h

#define MAX_STRING_LENGTH 1000000  // Increased for longer sequences
#define MAX_LINE_LENGTH 1024
#define MAX_SEQUENCES 1000

// Structure to store FASTA sequence
typedef struct {
    char* header;
    char* sequence;
} FastaSequence;

// Function to read FASTA file
int read_fasta(const char* filename, FastaSequence** sequences) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int seq_count = 0;
    char* current_sequence = NULL;
    int current_seq_len = 0;
    int current_seq_capacity = MAX_STRING_LENGTH;

    *sequences = (FastaSequence*)malloc(MAX_SEQUENCES * sizeof(FastaSequence));

    while (fgets(line, MAX_LINE_LENGTH, file)) {
        // Remove newline
        line[strcspn(line, "\n")] = 0;

        if (line[0] == '>') {
            // New sequence header
            if (current_sequence != NULL) {
                current_sequence[current_seq_len] = '\0';
                (*sequences)[seq_count - 1].sequence = current_sequence;
            }

            (*sequences)[seq_count].header = strdup(line + 1);
            current_sequence = (char*)malloc(current_seq_capacity * sizeof(char));
            current_seq_len = 0;
            seq_count++;
        } else {
            // Sequence content
            for (int i = 0; line[i]; i++) {
                if (!isspace(line[i])) {
                    current_sequence[current_seq_len++] = line[i];
                }
            }
        }
    }

    // Store last sequence
    if (current_sequence != NULL) {
        current_sequence[current_seq_len] = '\0';
        (*sequences)[seq_count - 1].sequence = current_sequence;
    }

    fclose(file);
    return seq_count;
}

int main(int argc, char* argv[]) {
    int rank, size;
    FastaSequence* sequences = NULL;
    int num_sequences;
    int k;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        // Master process reads input
        if (argc != 3) {
            printf("Usage: %s <fasta_file> <k>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        // Read FASTA file
        num_sequences = read_fasta(argv[1], &sequences);
        if (num_sequences < 0) {
            printf("Error reading FASTA file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        k = atoi(argv[2]);
        printf("Read %d sequences from FASTA file\n", num_sequences);
    }

    // Broadcast k and number of sequences to all processes
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_sequences, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // All processes allocate space for sequences
    if (rank != 0) {
        sequences = (FastaSequence*)malloc(num_sequences * sizeof(FastaSequence));
        for (int i = 0; i < num_sequences; i++) {
            sequences[i].sequence = (char*)malloc(MAX_STRING_LENGTH * sizeof(char));
            sequences[i].header = (char*)malloc(MAX_LINE_LENGTH * sizeof(char));
        }
    }

    // Broadcast sequences and headers to all processes
    for (int i = 0; i < num_sequences; i++) {
        MPI_Bcast(sequences[i].sequence, MAX_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(sequences[i].header, MAX_LINE_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    // Prepare to gather results
    char* local_result = (char*)malloc(100000 * sizeof(char));
    local_result[0] = '\0';  // Initialize as empty string
    char* global_result = NULL;

    if (rank == 0) {
        global_result = (char*)malloc(size * 100000 * sizeof(char));
    }

    // Distribute work among processes
    int total_pairs = (num_sequences * (num_sequences - 1)) / 2;
    for (int pair_idx = rank; pair_idx < total_pairs; pair_idx += size) {
        int i = 0;
        int j = 1;
        int remaining = pair_idx;
        while (remaining >= num_sequences - j) {
            remaining -= (num_sequences - j);
            i++;
            j = i + 1;
        }
        j += remaining;

        // Compute k-LCP for this pair
        int** LCP = compute_k_LCP(sequences[i].sequence, sequences[j].sequence, k);

        // Serialize result into local_result
        char buffer[1024];
        snprintf(buffer, sizeof(buffer), "LCP Table for %s and %s (k=%d):\n",
                 sequences[i].header, sequences[j].header, k);
        strcat(local_result, buffer);

        int n = strlen(sequences[i].sequence);
        int m = strlen(sequences[j].sequence);

        for (int row = 0; row < n; row++) {
            for (int col = 0; col < m; col++) {
                snprintf(buffer, sizeof(buffer), "%d ", LCP[row][col]);
                strcat(local_result, buffer);
            }
            strcat(local_result, "\n");
        }
        strcat(local_result, "\n");

        // Free the LCP table
        for (int row = 0; row < n; row++) {
            free(LCP[row]);
        }
        free(LCP);
    }

    // Gather all results to rank 0
    MPI_Gather(local_result, 100000, MPI_CHAR, global_result, 100000, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Rank 0 writes to file
    if (rank == 0) {
        FILE* result_file = fopen("result.txt", "w");
        if (!result_file) {
            printf("Error opening result file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        for (int i = 0; i < size; i++) {
            fprintf(result_file, "%s", &global_result[i * 100000]);
        }
        fclose(result_file);
        free(global_result);
    }

    // Clean up
    free(local_result);
    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i].sequence);
        free(sequences[i].header);
    }
    free(sequences);

    MPI_Finalize();
    return 0;
}