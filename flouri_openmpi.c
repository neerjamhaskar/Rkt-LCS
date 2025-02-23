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
int read_fake_fasta(const char* filename, FastaSequence** sequences) {
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

// Parallel version of Rkt_LCS using OpenMPI.
// Parallelization is done by distributing the work over the sequences S[0]..S[m-1].
int * Rkt_LCS_MPI(char* S[], int m, int k, int t) {
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Each process will compute a local best candidate.
    int local_len_max = 0;
    int local_start_index = -1;
    int local_pivot_index = -1;

    // Distribute the work over S[0]..S[m-1]:
    // Each process takes i values where i % size == rank.
    for (int i = rank; i < m; i += size) {
        int li = strlen(S[i]);
        const int numTables = m - i;
        int** LCP_i = (int**)malloc(numTables * sizeof(int*));
        if (!LCP_i) {
            fprintf(stderr, "Memory allocation failed for LCP_i\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int j = i; j < m; j++) {
            LCP_i[j - i] = compute_k_LCP_max(S[i], S[j], k);
        }
        
        // For each possible starting position p in S[i]:
        for (int p = 0; p < li; p++) {
            int** LengthStat = compute_LengthStat(LCP_i, S, m, p, i);
            int L = li - p;  // Maximum possible substring length from S[i] starting at p.
            
            // Check candidate lengths from longest to shortest.
            for (int l = L; l >= 1; l--) {  
                if (LengthStat[l - 1][m] >= t && l > local_len_max) {
                    local_start_index = p;
                    local_pivot_index = i;
                    local_len_max = l;
                }
            }
            free2DArray((int **)LengthStat, L);
        }
        free2DArray((int **)LCP_i, numTables);
    }
    
    // We use MPI_Allreduce to find the maximum candidate length overall.
    int global_len_max = 0;
    MPI_Allreduce(&local_len_max, &global_len_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                
    // Also gather the candidate start index.
    int gathered_starts[size];
    MPI_Gather(&local_start_index, 1, MPI_INT,
                gathered_starts, 1, MPI_INT,
                0, MPI_COMM_WORLD);

    // Also gather the pivot index.
    int gathered_pivot_index[size];
    MPI_Gather(&local_pivot_index, 1, MPI_INT,
                gathered_pivot_index, 1, MPI_INT,
                0, MPI_COMM_WORLD);
    
    // Also gather the candidate lengths.
    int gathered_lengths[size];
    MPI_Gather(&local_len_max, 1, MPI_INT,
               gathered_lengths, 1, MPI_INT,
               0, MPI_COMM_WORLD);
    
    // result: index i of pivot string, start index p in pivot string
    int* result = (int*)malloc(3 * sizeof(int));
    if (rank == 0) {
        // Rank 0 selects the candidate whose length equals global_len_max.
        // (If more than one candidate qualifies, we choose the first one.)
        for (int i = 0; i < size; i++) {
            if (gathered_lengths[i] == global_len_max && global_len_max > 0) {
                result[0] = gathered_pivot_index[i];
                result[1] = gathered_starts[i];
                break;
            }
        }
        // If no candidate was found, result remains NULL.
    }
    result[2] = global_len_max;
    
    return result;  // Only rank 0’s result is considered authoritative.
}

// Function to read a FASTQ-like file where only the (4k+2) header and (4k+3) sequence lines are used.
int read_fasta(const char* filename, FastaSequence** sequences) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int seq_count = 0;
    int line_num = 0;

    // Allocate array for sequences
    *sequences = (FastaSequence*)malloc(MAX_SEQUENCES * sizeof(FastaSequence));
    if (!*sequences) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(file);
        return -1;
    }

    while (fgets(line, MAX_LINE_LENGTH, file)) {
        line[strcspn(line, "\n")] = '\0';  // Remove newline

        if (line_num % 4 == 1) {
            // This is a header line (4k+2 when counting from 1)
            (*sequences)[seq_count].header = strdup(line);
            if (!(*sequences)[seq_count].header) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                fclose(file);
                return -1;
            }
        } else if (line_num % 4 == 2) {
            // This is the sequence line (4k+3 when counting from 1)
            (*sequences)[seq_count].sequence = strdup(line);
            if (!(*sequences)[seq_count].sequence) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                fclose(file);
                return -1;
            }
            seq_count++; // Increment count after both header and sequence have been read.
        }
        // Lines where (line_num % 4 == 0) or (line_num % 4 == 3) are ignored.
        line_num++;
    }

    fclose(file);
    return seq_count;
}

// Main function that reads a FASTA file and calls Rkt_LCS_MPI.
int main(int argc, char* argv[]) {
    int rank, size;
    FastaSequence* sequences = NULL;
    int num_sequences;
    int k, t;  // k, t: args for Rkt_LCS_MPI

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Rank 0 reads the FASTA file and parses command-line arguments.
    if (rank == 0) {
        if (argc != 4) {
            printf("Usage: %s <fasta_file> <k> <t>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        num_sequences = read_fasta(argv[1], &sequences);
        if (num_sequences < 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        k = atoi(argv[2]);
        t = atoi(argv[3]);
        printf("Read %d sequences from FASTA file\n", num_sequences);
    }

    // Broadcast k, t, and the number of sequences to all processes.
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_sequences, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // All processes allocate memory for the sequences if they are not rank 0.
    if (rank != 0) {
        sequences = (FastaSequence*)malloc(num_sequences * sizeof(FastaSequence));
        if (!sequences) {
            fprintf(stderr, "Memory allocation failed for sequences on rank %d\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < num_sequences; i++) {
            sequences[i].header = (char*)malloc(MAX_LINE_LENGTH);
            sequences[i].sequence = (char*)malloc(MAX_STRING_LENGTH);
            if (!sequences[i].header || !sequences[i].sequence) {
                fprintf(stderr, "Memory allocation failed for sequence %d on rank %d\n", i, rank);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    // Broadcast each sequence's header and sequence from rank 0 to all processes.
    for (int i = 0; i < num_sequences; i++) {
        MPI_Bcast(sequences[i].header, MAX_LINE_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(sequences[i].sequence, MAX_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    // Create an array of char* pointers (one per sequence) that points to each sequence string.
    char** S_array = (char**)malloc(num_sequences * sizeof(char*));
    if (!S_array) {
        fprintf(stderr, "Memory allocation failed for S_array on rank %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < num_sequences; i++) {
        S_array[i] = sequences[i].sequence;
    }

    // Call the MPI-parallelized longest common substring function.
    // This function will run in parallel over all processes and return the best candidate.
    int* lcs_result = Rkt_LCS_MPI(S_array, num_sequences, k, t);

    // Rank 0 prints the result.
    if (rank == 0) {
        if (lcs_result != NULL && lcs_result[2] > 0) {
            int pivot_index = lcs_result[0];
            int start_index = lcs_result[1];
            int length = lcs_result[2];
            printf("Best Rk-tLCS is in string %d, starting at index %d: %s\n", 
                lcs_result[0], 
                lcs_result[1], 
                substring(S_array[pivot_index], start_index, length));
            free(lcs_result);
        } else {
            printf("No valid LCS found.\n");
        }
    }

    // Free allocated memory.
    free(S_array);
    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i].sequence);
        free(sequences[i].header);
    }
    free(sequences);

    MPI_Finalize();
    return 0;
}




/*** main function that calls compute_k_LCP with openmpi ***/
// int main(int argc, char* argv[]) {
//     int rank, size;
//     FastaSequence* sequences = NULL;
//     int num_sequences;
//     int k;

//     MPI_Init(&argc, &argv);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     if (rank == 0) {
//         if (argc != 3) {
//             printf("Usage: %s <fasta_file> <k>\n", argv[0]);
//             MPI_Abort(MPI_COMM_WORLD, 1);
//         }

//         num_sequences = read_fasta(argv[1], &sequences);
//         if (num_sequences < 0) {
//             MPI_Abort(MPI_COMM_WORLD, 1);
//         }

//         k = atoi(argv[2]);
//         printf("Read %d sequences from FASTA file\n", num_sequences);
//     }

//     MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(&num_sequences, 1, MPI_INT, 0, MPI_COMM_WORLD);

//     if (rank != 0) {
//         sequences = (FastaSequence*)malloc(num_sequences * sizeof(FastaSequence));
//         for (int i = 0; i < num_sequences; i++) {
//             sequences[i].header = (char*)malloc(MAX_LINE_LENGTH);
//             sequences[i].sequence = (char*)malloc(MAX_STRING_LENGTH);
//         }
//     }

//     for (int i = 0; i < num_sequences; i++) {
//         MPI_Bcast(sequences[i].header, MAX_LINE_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
//         MPI_Bcast(sequences[i].sequence, MAX_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
//     }

//     char* local_result = (char*)malloc(BUFFER_SIZE);
//     local_result[0] = '\0';
//     char* global_result = NULL;

//     if (rank == 0) {
//         global_result = (char*)malloc(size * BUFFER_SIZE);
//     }

//     int total_pairs = (num_sequences * (num_sequences - 1)) / 2;
//     for (int pair_idx = rank; pair_idx < total_pairs; pair_idx += size) {
//         int i = 0, j = 1, remaining = pair_idx;
//         while (remaining >= num_sequences - j) {
//             remaining -= (num_sequences - j);
//             i++;
//             j = i + 1;
//         }
//         j += remaining;

//         int** LCP = compute_k_LCP(sequences[i].sequence, sequences[j].sequence, k);
//         char buffer[BUFFER_SIZE];
//         snprintf(buffer, BUFFER_SIZE, "LCP Table for %s and %s:\n", sequences[i].header, sequences[j].header);
//         strncat(local_result, buffer, BUFFER_SIZE);

//         for (int row = 0; row < strlen(sequences[i].sequence); row++) {
//             for (int col = 0; col < strlen(sequences[j].sequence); col++) {
//                 snprintf(buffer, BUFFER_SIZE, "%d ", LCP[row][col]);
//                 strncat(local_result, buffer, BUFFER_SIZE);
//             }
//             strncat(local_result, "\n", BUFFER_SIZE);
//         }

//         for (int row = 0; row < strlen(sequences[i].sequence); row++) {
//             free(LCP[row]);
//         }
//         free(LCP);
//     }

//     MPI_Gather(local_result, BUFFER_SIZE, MPI_CHAR, global_result, BUFFER_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

//     if (rank == 0) {
//         FILE* result_file = fopen("result.txt", "w");
//         if (!result_file) {
//             printf("Error opening result file\n");
//             MPI_Abort(MPI_COMM_WORLD, 1);
//         }

//         for (int i = 0; i < size; i++) {
//             fprintf(result_file, "%s", &global_result[i * BUFFER_SIZE]);
//         }
//         fclose(result_file);
//         free(global_result);
//     }

//     free(local_result);
//     for (int i = 0; i < num_sequences; i++) {
//         free(sequences[i].sequence);
//         free(sequences[i].header);
//     }
//     free(sequences);

//     MPI_Finalize();
//     return 0;
// }