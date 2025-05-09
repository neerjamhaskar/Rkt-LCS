#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "flouri_cpu.h"

// Structure to hold the result for each sequence
typedef struct {
    int seq_num;          // Original sequence number (1-based)
    PosLenKey result;     // Result from Rkt_LCS_single
    char substring[1024]; // The actual substring found
} SequenceResult;

int main(int argc, char* argv[]) {
    int rank, size;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check command line arguments (only in rank 0)
    if (rank == 0) {
        printf("DEBUG: Starting program with %d processes\n", size);
        if (argc != 5) {
            printf("Usage: %s <fasta_file> <k> <t> <tau>\n", argv[0]);
            printf("  fasta_file: Input FASTA file containing sequences\n");
            printf("  k: Maximum number of mismatches allowed\n");
            printf("  t: Minimum number of S2 strings that must contain the substring\n");
            printf("  tau: Minimum match length\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
    }

    // Parse command line arguments
    const char* fasta_file = argv[1];
    int k = atoi(argv[2]);
    int t = atoi(argv[3]);
    int tau = atoi(argv[4]);

    if (rank == 0) {
        printf("DEBUG: Command arguments - file: %s, k: %d, t: %d, tau: %d\n", 
               fasta_file, k, t, tau);
    }

    // Read sequences from FASTA file (all ranks need this)
    int num_sequences = 0;
    if (rank == 0) printf("DEBUG: Reading FASTA file %s\n", fasta_file);
    
    FastaSequence* sequences = read_fasta(fasta_file, &num_sequences);
    
    if (rank == 0) printf("DEBUG: Read %d sequences\n", num_sequences);
    
    if (!sequences) {
        if (rank == 0) printf("Error: Failed to read sequences from %s\n", fasta_file);
        MPI_Finalize();
        return 1;
    }

    // Basic validation
    if (num_sequences < 2 || k < 0 || t < 1 || tau < 1) {
        if (rank == 0) {
            printf("Error: Invalid parameters or insufficient sequences\n");
        }
        free_fasta(sequences, num_sequences);
        MPI_Finalize();
        return 1;
    }

    // Calculate work distribution
    int sequences_per_proc = num_sequences / size;
    int extra_sequences = num_sequences % size;
    int start_seq = rank * sequences_per_proc + (rank < extra_sequences ? rank : extra_sequences);
    int end_seq = start_seq + sequences_per_proc + (rank < extra_sequences ? 1 : 0);

    if (rank == 0) {
        printf("DEBUG: Work distribution - total sequences: %d\n", num_sequences);
        for (int i = 0; i < size; i++) {
            int start = i * sequences_per_proc + (i < extra_sequences ? i : extra_sequences);
            int end = start + sequences_per_proc + (i < extra_sequences ? 1 : 0);
            printf("DEBUG: Process %d will handle sequences %d to %d\n", i, start, end-1);
        }
    }

    // Allocate array for local results
    int local_count = end_seq - start_seq;
    printf("Process %d: Processing %d sequences (%d to %d)\n", 
           rank, local_count, start_seq, end_seq-1);
    
    SequenceResult* local_results = malloc(local_count * sizeof(SequenceResult));
    if (!local_results) {
        printf("Error: Memory allocation failed for local_results on rank %d\n", rank);
        free_fasta(sequences, num_sequences);
        MPI_Finalize();
        return 1;
    }
    
    // Process assigned sequences
    for (int i = 0; i < local_count; i++) {
        int seq_idx = start_seq + i;
        // printf("Process %d: Processing sequence %d (%d of %d)\n", 
        //        rank, seq_idx, i+1, local_count);
               
        const char* S1 = sequences[seq_idx].sequence;
        int r = strlen(S1);
        
        // printf("Process %d: S1 length = %d\n", rank, r);

        // Create S2_array with all sequences except S1
        char** S2_array = malloc((num_sequences - 1) * sizeof(char*));
        if (!S2_array) {
            printf("Error: Memory allocation failed for S2_array on rank %d\n", rank);
            free(local_results);
            free_fasta(sequences, num_sequences);
            MPI_Finalize();
            return 1;
        }
        
        int s2_idx = 0;
        for (int j = 0; j < num_sequences; j++) {
            if (j != seq_idx) {
                S2_array[s2_idx++] = sequences[j].sequence;
            }
        }

        // printf("Process %d: Calling Rkt_LCS_single (k=%d, t=%d, tau=%d, r=%d)\n", 
        //        rank, k, t, tau, r);
               
        // Compute result for this sequence
        PosLenKey result = Rkt_LCS_single(S1, S2_array, num_sequences - 1, k, t, tau, r);

        // printf("Process %d: Got result position=%d, length=%d\n", 
        //        rank, result.position, result.length);
               
        // Store result
        local_results[i].seq_num = seq_idx + 1;  // 1-based sequence number
        local_results[i].result = result;
        
        // Store the substring if a valid result was found
        if (result.position != -1 && result.length > 0) {
            strncpy(local_results[i].substring, S1 + result.position, result.length);
            local_results[i].substring[result.length] = '\0';
            // printf("Process %d: Substring = %s\n", rank, local_results[i].substring);
        } else {
            local_results[i].substring[0] = '\0';
            // printf("Process %d: No valid substring found\n", rank);
        }

        free(S2_array);
    }

    // printf("Process %d: All local processing complete\n", rank);

    // Gather the number of results from each process
    int* result_counts = NULL;
    if (rank == 0) {
        result_counts = malloc(size * sizeof(int));
    }
    
    MPI_Gather(&local_count, 1, MPI_INT, result_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Calculate displacements for Gatherv
    int* displacements = NULL;
    int total_results = 0;
    if (rank == 0) {
        displacements = malloc(size * sizeof(int));
        displacements[0] = 0;
        
        for (int i = 0; i < size; i++) {
            if (i > 0) {
                displacements[i] = displacements[i-1] + result_counts[i-1];
            }
            total_results += result_counts[i];
        }
        
        printf("Total results across all processes: %d\n", total_results);
    }
    
    // Create MPI datatype for SequenceResult
    MPI_Datatype result_type;
    int blocklengths[3] = {1, 2, 1024};
    MPI_Aint offsets[3];
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_CHAR};
    
    offsets[0] = offsetof(SequenceResult, seq_num);
    offsets[1] = offsetof(SequenceResult, result);
    offsets[2] = offsetof(SequenceResult, substring);
    
    MPI_Type_create_struct(3, blocklengths, offsets, types, &result_type);
    MPI_Type_commit(&result_type);
    
    // Gather all results to rank 0
    SequenceResult* all_results = NULL;
    if (rank == 0) {
        all_results = malloc(total_results * sizeof(SequenceResult));
    }
    
    MPI_Gatherv(local_results, local_count, result_type, 
                all_results, result_counts, displacements, 
                result_type, 0, MPI_COMM_WORLD);
    
    // Process 0 writes the combined results
    if (rank == 0) {
        FILE* outfile = fopen("rkt_lcs_cpu_results.txt", "w");
        if (outfile) {
            fprintf(outfile, "Combined results with k=%d, t=%d, tau=%d:\n", k, t, tau);
            fprintf(outfile, "===============================================\n\n");
            
            for (int i = 0; i < total_results; i++) {
                int seq_num = all_results[i].seq_num;
                fprintf(outfile, "Sequence %d (%s):\n", 
                        seq_num, 
                        sequences[seq_num - 1].name);
                
                if (all_results[i].result.position != -1) {
                    fprintf(outfile, "  Position: %d\n", all_results[i].result.position);
                    fprintf(outfile, "  Length: %d\n", all_results[i].result.length);
                    fprintf(outfile, "  Substring: %s\n", all_results[i].substring);
                } else {
                    fprintf(outfile, "  No valid substring found\n");
                }
                fprintf(outfile, "----------------------------------------\n");
            }
            
            fclose(outfile);
            printf("All results written to rkt_lcs_results.txt\n");
            fflush(stdout);
            
            free(all_results);
            free(result_counts);
            free(displacements);
        } else {
            printf("Error: Could not create results file\n");
        }
    }
    
    // Clean up
    free(local_results);
    free_fasta(sequences, num_sequences);
    MPI_Type_free(&result_type);

    
    // printf("Process %d: Finalizing\n", rank);
    MPI_Finalize();
    return 0;
}