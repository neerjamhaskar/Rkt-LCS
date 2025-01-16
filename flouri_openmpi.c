////  ./dependencies/bin/mpicc -o flouri_openmpi.o flouri_openmpi.c
///// ./dependencies/bin/mpirun -np 4 ./flouri_openmpi.o fake_reads.fasta 3
#include "dependencies/include/mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "flouri_lcp_table.h" // Include table.h

#define MAX_STRING_LENGTH 1000000 // Maximum sequence length
#define MAX_LINE_LENGTH 1024
#define MAX_SEQUENCES 1000
#define BUFFER_SIZE 1024 // To avoid strcat overflows

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
    if (!*sequences) {
        printf("Error: Memory allocation failed\n");
        fclose(file);
        return -1;
    }

    while (fgets(line, MAX_LINE_LENGTH, file)) {
        line[strcspn(line, "\n")] = 0; // Remove newline

        if (line[0] == '>') {
            if (current_sequence) {
                current_sequence[current_seq_len] = '\0';
                (*sequences)[seq_count - 1].sequence = current_sequence;
            }

            (*sequences)[seq_count].header = strdup(line + 1);
            current_sequence = (char*)malloc(current_seq_capacity * sizeof(char));
            if (!current_sequence || !(*sequences)[seq_count].header) {
                printf("Error: Memory allocation failed\n");
                fclose(file);
                return -1;
            }

            current_seq_len = 0;
            seq_count++;
        } else {
            for (int i = 0; line[i]; i++) {
                if (!isspace(line[i])) {
                    current_sequence[current_seq_len++] = line[i];
                }
            }
        }
    }

    if (current_sequence) {
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
        if (argc != 3) {
            printf("Usage: %s <fasta_file> <k>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        num_sequences = read_fasta(argv[1], &sequences);
        if (num_sequences < 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        k = atoi(argv[2]);
        printf("Read %d sequences from FASTA file\n", num_sequences);
    }

    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_sequences, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        sequences = (FastaSequence*)malloc(num_sequences * sizeof(FastaSequence));
        for (int i = 0; i < num_sequences; i++) {
            sequences[i].header = (char*)malloc(MAX_LINE_LENGTH);
            sequences[i].sequence = (char*)malloc(MAX_STRING_LENGTH);
        }
    }

    for (int i = 0; i < num_sequences; i++) {
        MPI_Bcast(sequences[i].header, MAX_LINE_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(sequences[i].sequence, MAX_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    char* local_result = (char*)malloc(BUFFER_SIZE);
    local_result[0] = '\0';
    char* global_result = NULL;

    if (rank == 0) {
        global_result = (char*)malloc(size * BUFFER_SIZE);
    }

    int total_pairs = (num_sequences * (num_sequences - 1)) / 2;
    for (int pair_idx = rank; pair_idx < total_pairs; pair_idx += size) {
        int i = 0, j = 1, remaining = pair_idx;
        while (remaining >= num_sequences - j) {
            remaining -= (num_sequences - j);
            i++;
            j = i + 1;
        }
        j += remaining;

        int** LCP = compute_k_LCP(sequences[i].sequence, sequences[j].sequence, k);
        char buffer[BUFFER_SIZE];
        snprintf(buffer, BUFFER_SIZE, "LCP Table for %s and %s:\n", sequences[i].header, sequences[j].header);
        strncat(local_result, buffer, BUFFER_SIZE);

        for (int row = 0; row < strlen(sequences[i].sequence); row++) {
            for (int col = 0; col < strlen(sequences[j].sequence); col++) {
                snprintf(buffer, BUFFER_SIZE, "%d ", LCP[row][col]);
                strncat(local_result, buffer, BUFFER_SIZE);
            }
            strncat(local_result, "\n", BUFFER_SIZE);
        }

        for (int row = 0; row < strlen(sequences[i].sequence); row++) {
            free(LCP[row]);
        }
        free(LCP);
    }

    MPI_Gather(local_result, BUFFER_SIZE, MPI_CHAR, global_result, BUFFER_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        FILE* result_file = fopen("result.txt", "w");
        if (!result_file) {
            printf("Error opening result file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < size; i++) {
            fprintf(result_file, "%s", &global_result[i * BUFFER_SIZE]);
        }
        fclose(result_file);
        free(global_result);
    }

    free(local_result);
    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i].sequence);
        free(sequences[i].header);
    }
    free(sequences);

    MPI_Finalize();
    return 0;
}