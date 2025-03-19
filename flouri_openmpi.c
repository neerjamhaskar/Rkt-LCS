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
#define MAX_SEQUENCES 10000 // Increased to handle larger datasets
#define BUFFER_SIZE (1024 * 1024) // 1MB buffer
#define MAX_RESULT_SIZE (BUFFER_SIZE * 32) // Maximum size for global result
#define CHUNK_SIZE 50 // Reduced chunk size for better work distribution
#define MAX_PAIRS_PER_CHUNK 1000 // Maximum number of pairs to process in each batch

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
        printf("Total pairs to process: %lld\n", (long long)num_sequences * (num_sequences - 1) / 2);
    }

    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_sequences, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Process sequences in chunks to manage memory better
    int num_chunks = (num_sequences + CHUNK_SIZE - 1) / CHUNK_SIZE;
    printf("Rank %d: Processing %d chunks\n", rank, num_chunks);

    // Allocate memory for two chunks (current and target)
    FastaSequence* current_chunk = (FastaSequence*)malloc(CHUNK_SIZE * sizeof(FastaSequence));
    FastaSequence* target_chunk = (FastaSequence*)malloc(CHUNK_SIZE * sizeof(FastaSequence));
    
    if (!current_chunk || !target_chunk) {
        printf("Rank %d: Failed to allocate chunk arrays\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize chunks
    for (int i = 0; i < CHUNK_SIZE; i++) {
        current_chunk[i].header = (char*)malloc(MAX_LINE_LENGTH);
        current_chunk[i].sequence = (char*)malloc(MAX_STRING_LENGTH);
        target_chunk[i].header = (char*)malloc(MAX_LINE_LENGTH);
        target_chunk[i].sequence = (char*)malloc(MAX_STRING_LENGTH);
        
        if (!current_chunk[i].header || !current_chunk[i].sequence ||
            !target_chunk[i].header || !target_chunk[i].sequence) {
            printf("Rank %d: Failed to allocate sequence memory\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Allocate and initialize local result buffer
    char* local_result = (char*)calloc(BUFFER_SIZE, sizeof(char));
    if (!local_result) {
        printf("Rank %d: Failed to allocate local_result\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    local_result[0] = '\0';

    char* global_result = NULL;
    if (rank == 0) {
        global_result = (char*)calloc(MAX_RESULT_SIZE, sizeof(char));
        if (!global_result) {
            printf("Rank 0: Failed to allocate global_result\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    size_t current_pos = 0;
    long long total_pairs_processed = 0;

    // Process all pairs between chunks
    for (int i_chunk = 0; i_chunk < num_chunks; i_chunk++) {
        int i_start = i_chunk * CHUNK_SIZE;
        int i_end = (i_chunk + 1) * CHUNK_SIZE;
        if (i_end > num_sequences) i_end = num_sequences;
        int i_chunk_size = i_end - i_start;

        if (rank == 0) {
            printf("Processing chunk %d/%d (sequences %d to %d)\n", 
                   i_chunk + 1, num_chunks, i_start, i_end - 1);
        }

        // Load current chunk
        if (rank == 0) {
            for (int i = 0; i < i_chunk_size; i++) {
                strcpy(current_chunk[i].header, sequences[i_start + i].header);
                strcpy(current_chunk[i].sequence, sequences[i_start + i].sequence);
            }
        }

        // Broadcast current chunk to all processes
        for (int i = 0; i < i_chunk_size; i++) {
            MPI_Bcast(current_chunk[i].header, MAX_LINE_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Bcast(current_chunk[i].sequence, MAX_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        }

        // Process pairs with all subsequent chunks
        for (int j_chunk = i_chunk; j_chunk < num_chunks; j_chunk++) {
            int j_start = j_chunk * CHUNK_SIZE;
            int j_end = (j_chunk + 1) * CHUNK_SIZE;
            if (j_end > num_sequences) j_end = num_sequences;
            int j_chunk_size = j_end - j_start;

            // Load target chunk
            if (rank == 0) {
                for (int j = 0; j < j_chunk_size; j++) {
                    strcpy(target_chunk[j].header, sequences[j_start + j].header);
                    strcpy(target_chunk[j].sequence, sequences[j_start + j].sequence);
                }
            }

            // Broadcast target chunk to all processes
            for (int j = 0; j < j_chunk_size; j++) {
                MPI_Bcast(target_chunk[j].header, MAX_LINE_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
                MPI_Bcast(target_chunk[j].sequence, MAX_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
            }

            // Calculate pairs for this chunk combination
            int total_pairs;
            if (i_chunk == j_chunk) {
                // Within same chunk: compute pairs (n choose 2)
                total_pairs = (i_chunk_size * (i_chunk_size - 1)) / 2;
            } else {
                // Between different chunks: compute all pairs
                total_pairs = i_chunk_size * j_chunk_size;
            }

            // Distribute work among processes
            int pairs_per_process = (total_pairs + size - 1) / size;
            int my_start = rank * pairs_per_process;
            int my_end = (rank + 1) * pairs_per_process;
            if (my_end > total_pairs) my_end = total_pairs;

            if (rank == 0) {
                printf("Processing %d pairs between chunks %d and %d\n", total_pairs, i_chunk, j_chunk);
            }

            // Process pairs in smaller batches
            for (int batch_start = my_start; batch_start < my_end; batch_start += MAX_PAIRS_PER_CHUNK) {
                int batch_end = batch_start + MAX_PAIRS_PER_CHUNK;
                if (batch_end > my_end) batch_end = my_end;

                for (int pair_idx = batch_start; pair_idx < batch_end; pair_idx++) {
                    int i, j;
                    if (i_chunk == j_chunk) {
                        // Within same chunk: use triangular number formula
                        i = 0;
                        j = 1;
                        int remaining = pair_idx;
                        while (remaining >= i_chunk_size - j) {
                            remaining -= (i_chunk_size - j);
                            i++;
                            j = i + 1;
                        }
                        j += remaining;
                    } else {
                        // Between different chunks: use division and modulo
                        i = pair_idx / j_chunk_size;
                        j = pair_idx % j_chunk_size;
                    }

                    int** LCP = compute_k_LCP(current_chunk[i].sequence, target_chunk[j].sequence, k);
                    if (!LCP) {
                        printf("Rank %d: Failed to compute LCP for pair %d-%d\n", rank, i_start + i, j_start + j);
                        continue;
                    }

                    char buffer[BUFFER_SIZE];
                    int written = snprintf(buffer, BUFFER_SIZE, "LCP Table for %s and %s:\n", 
                                         current_chunk[i].header, target_chunk[j].header);
                    if (written > 0 && current_pos + written < BUFFER_SIZE) {
                        strcat(local_result + current_pos, buffer);
                        current_pos += written;
                    }

                    for (int row = 0; row < strlen(current_chunk[i].sequence); row++) {
                        for (int col = 0; col < strlen(target_chunk[j].sequence); col++) {
                            written = snprintf(buffer, BUFFER_SIZE, "%d ", LCP[row][col]);
                            if (written > 0 && current_pos + written < BUFFER_SIZE) {
                                strcat(local_result + current_pos, buffer);
                                current_pos += written;
                            }
                        }
                        if (current_pos + 1 < BUFFER_SIZE) {
                            strcat(local_result + current_pos, "\n");
                            current_pos++;
                        }
                    }

                    // Free LCP table
                    for (int row = 0; row < strlen(current_chunk[i].sequence); row++) {
                        free(LCP[row]);
                    }
                    free(LCP);

                    total_pairs_processed++;
                    if (rank == 0 && total_pairs_processed % 1000 == 0) {
                        printf("Processed %lld pairs\n", total_pairs_processed);
                    }
                }

                // Gather results for this batch
                MPI_Gather(local_result, BUFFER_SIZE, MPI_CHAR, 
                          global_result, BUFFER_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

                if (rank == 0) {
                    FILE* result_file = fopen("result.txt", "a");
                    if (!result_file) {
                        printf("Error opening result file\n");
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }

                    fprintf(result_file, "%s", global_result);
                    fclose(result_file);
                }

                // Reset local result buffer for next batch
                memset(local_result, 0, BUFFER_SIZE);
                current_pos = 0;
            }
        }
    }

    // Sum up total pairs processed across all processes
    long long global_total_pairs;
    MPI_Reduce(&total_pairs_processed, &global_total_pairs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Total pairs processed: %lld\n", global_total_pairs);
    }

    // Cleanup
    free(local_result);
    if (rank == 0) {
        free(global_result);
    }
    for (int i = 0; i < CHUNK_SIZE; i++) {
        free(current_chunk[i].sequence);
        free(current_chunk[i].header);
        free(target_chunk[i].sequence);
        free(target_chunk[i].header);
    }
    free(current_chunk);
    free(target_chunk);
    if (rank == 0) {
        for (int i = 0; i < num_sequences; i++) {
            free(sequences[i].sequence);
            free(sequences[i].header);
        }
        free(sequences);
    }

    MPI_Finalize();
    return 0;
}