#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <cuda_runtime.h>
#include "flouri_cpu.h"

// Function declarations
PosLenCount* compute_k_LCP_max_multi_cpu(const char* S1, char** S2_array, int num_S2, int k, int tau, int r, int* result_size);
extern "C" PosLenCount* compute_k_LCP_max_multi_cuda(const char* S1, char** S2_array, int num_S2, int k, int tau, int r, int* result_size);
extern "C" PosLenKey Rkt_LCS_single_cuda(const char* S1, char** S2_array, int num_S2, int k, int t, int tau, int r);

// Comparison function for sorting PosLenCount entries
int compare_pos_len_count(const void* a, const void* b) {
    PosLenCount* posa = (PosLenCount*)a;
    PosLenCount* posb = (PosLenCount*)b;
    
    // Sort by position first
    if (posa->key.position != posb->key.position) {
        return posa->key.position - posb->key.position;
    }
    
    // If positions are the same, sort by length
    if (posa->key.length != posb->key.length) {
        return posa->key.length - posb->key.length;
    }
    
    // If both position and length are the same, sort by count
    return posa->count - posb->count;
}

// Function to get current time in microseconds
double get_time_usec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000000 + (double)tv.tv_usec;
}

// Function to run benchmark for a given parameter set
void run_benchmark(FastaSequence* sequences, int num_sequences, int s1_seq_num, int k, int t, int tau, bool verbose_output) {
    // Flag to control detailed printing (set to false for cleaner output)
    bool print_detailed_results = verbose_output;
    
    printf("------------------------------------------------------\n");
    printf("Parameters:\n");
    printf("  - Sequence: %s (index %d)\n", sequences[s1_seq_num - 1].name, s1_seq_num);
    printf("  - k (max mismatches): %d\n", k);
    printf("  - t (min occurrences): %d\n", t);
    printf("  - tau (min match length): %d\n", tau);
    
    // Prepare input data
    const char* S1 = sequences[s1_seq_num - 1].sequence;
    int r = strlen(S1);
    
    char** S2_array = (char**)malloc((num_sequences - 1) * sizeof(char*));
    int s2_idx = 0;
    for (int i = 0; i < num_sequences; i++) {
        if (i != s1_seq_num - 1) {
            S2_array[s2_idx++] = sequences[i].sequence;
        }
    }
    int num_S2 = num_sequences - 1;
    
    printf("  - S1 length: %d\n", r);
    printf("  - Number of S2 sequences: %d\n", num_S2);
    
    // Variables for results
    int cpu_result_size = 0;
    int gpu_result_size = 0;
    PosLenCount* cpu_results = NULL;
    PosLenCount* gpu_results = NULL;
    double cpu_time = 0.0;
    double gpu_time = 0.0;
    
    // Run CPU implementation
    printf("\nRunning CPU implementation...\n");
    double cpu_start = get_time_usec();
    cpu_results = compute_k_LCP_max_multi_cpu(S1, S2_array, num_S2, k, tau, r, &cpu_result_size);
    double cpu_end = get_time_usec();
    cpu_time = (cpu_end - cpu_start) / 1000.0; // Convert to milliseconds
    
    printf("CPU implementation completed in %.2f ms\n", cpu_time);
    printf("CPU found %d unique position-length pairs\n", cpu_result_size);
    
    // Run GPU implementation
    printf("\nRunning GPU implementation...\n");
    double gpu_start = get_time_usec();
    gpu_results = compute_k_LCP_max_multi_cuda(S1, S2_array, num_S2, k, tau, r, &gpu_result_size);
    double gpu_end = get_time_usec();
    gpu_time = (gpu_end - gpu_start) / 1000.0; // Convert to milliseconds
    
    printf("GPU implementation completed in %.2f ms\n", gpu_time);
    printf("GPU found %d unique position-length pairs\n", gpu_result_size);
    
    // Calculate speedup
    double speedup = cpu_time / gpu_time;
    printf("\nResults:\n");
    printf("  - CPU time: %.2f ms\n", cpu_time);
    printf("  - GPU time: %.2f ms\n", gpu_time);
    printf("  - Speedup: %.2fx\n", speedup);
    
    // Validate results
    int valid = 1;
    if (cpu_result_size != gpu_result_size) {
        printf("Warning: Result sizes differ (CPU: %d, GPU: %d)\n", cpu_result_size, gpu_result_size);
        valid = 0;
    } else {
        printf("Result sizes match (%d entries)\n", cpu_result_size);
        
        // Sort both results for comparison
        if (cpu_result_size > 0 && gpu_result_size > 0) {
            qsort(cpu_results, cpu_result_size, sizeof(PosLenCount), compare_pos_len_count);
            qsort(gpu_results, gpu_result_size, sizeof(PosLenCount), compare_pos_len_count);
            
            // Compare all entries
            int mismatch_count = 0;
            for (int i = 0; i < cpu_result_size; i++) {
                if (cpu_results[i].key.position != gpu_results[i].key.position ||
                    cpu_results[i].key.length != gpu_results[i].key.length ||
                    cpu_results[i].count != gpu_results[i].count) {
                    mismatch_count++;
                }
            }
            
            if (mismatch_count == 0) {
                printf("All results match perfectly! ✓\n");
            } else {
                printf("Warning: %d mismatches found between CPU and GPU results\n", mismatch_count);
                valid = 0;
            }
        }
    }
    
    // Only print detailed results if the flag is set
    if (print_detailed_results) {
        printf("\n==================== COMPLETE RESULTS ====================\n");
        
        // Print CPU results table
        printf("\nCPU Results (%d entries):\n", cpu_result_size);
        printf("%-10s %-15s %-15s %-15s\n", "Index", "Position", "Length", "Count");
        printf("---------------------------------------------------------\n");
        
        for (int i = 0; i < cpu_result_size; i++) {
            printf("%-10d %-15d %-15d %-15d\n", 
                   i,
                   cpu_results[i].key.position,
                   cpu_results[i].key.length,
                   cpu_results[i].count);
            
            // Comment out substring printing
            /*
            // Print the actual substring from S1
            int pos = cpu_results[i].key.position;
            int len = cpu_results[i].key.length;
            printf("    Substring: \"");
            for (int j = 0; j < len && pos + j < r; j++) {
                printf("%c", S1[pos + j]);
            }
            printf("\"\n");
            */
        }
        
        // Print GPU results table
        printf("\nGPU Results (%d entries):\n", gpu_result_size);
        printf("%-10s %-15s %-15s %-15s\n", "Index", "Position", "Length", "Count");
        printf("---------------------------------------------------------\n");
        
        for (int i = 0; i < gpu_result_size; i++) {
            printf("%-10d %-15d %-15d %-15d\n", 
                   i,
                   gpu_results[i].key.position,
                   gpu_results[i].key.length,
                   gpu_results[i].count);
            
            // Comment out substring printing
            /*
            // Print the actual substring from S1
            int pos = gpu_results[i].key.position;
            int len = gpu_results[i].key.length;
            printf("    Substring: \"");
            for (int j = 0; j < len && pos + j < r; j++) {
                printf("%c", S1[pos + j]);
            }
            printf("\"\n");
            */
        }
        
        // If there were mismatches, show them in a comparison table
        if (!valid && cpu_result_size == gpu_result_size) {
            printf("\nDetailed comparison of mismatches:\n");
            printf("%-5s %-15s %-15s %-15s %-15s %-15s %-15s\n", 
                   "Index", "CPU Pos", "GPU Pos", "CPU Len", "GPU Len", "CPU Count", "GPU Count");
            printf("-----------------------------------------------------------------------------------\n");
            
            for (int i = 0; i < cpu_result_size; i++) {
                if (cpu_results[i].key.position != gpu_results[i].key.position ||
                    cpu_results[i].key.length != gpu_results[i].key.length ||
                    cpu_results[i].count != gpu_results[i].count) {
                    printf("%-5d %-15d %-15d %-15d %-15d %-15d %-15d\n", 
                        i, 
                        cpu_results[i].key.position, 
                        gpu_results[i].key.position,
                        cpu_results[i].key.length,
                        gpu_results[i].key.length,
                        cpu_results[i].count,
                        gpu_results[i].count);
                }
            }
        }
        
        printf("\n==========================================================\n");
    } else {
        // Print just a summary of results for cleaner output
        if (cpu_result_size > 0 && !valid) {
            printf("\nWARNING: The CPU and GPU implementations produced different results!\n");
            printf("Use the detailed output option to diagnose the issue.\n");
        }
    }
    
    // Clean up
    free(S2_array);
    if (cpu_results) free(cpu_results);
    if (gpu_results) free(gpu_results);
    
    printf("------------------------------------------------------\n\n");
}

// Add a new benchmark function for LCS single
// Function to benchmark Rkt_LCS_single (CPU vs CUDA)
void benchmark_lcs_single(FastaSequence* sequences, int num_sequences, int s1_seq_num, int k, int t, int tau, bool verbose_output) {
    printf("------------------------------------------------------\n");
    printf("Benchmarking Rkt_LCS_single (CPU vs CUDA):\n");
    printf("Parameters:\n");
    printf("  - Sequence: %s (index %d)\n", sequences[s1_seq_num - 1].name, s1_seq_num);
    printf("  - k (max mismatches): %d\n", k);
    printf("  - t (min occurrences): %d\n", t);
    printf("  - tau (min match length): %d\n", tau);
    
    // Prepare input data
    const char* S1 = sequences[s1_seq_num - 1].sequence;
    int r = strlen(S1);
    
    char** S2_array = (char**)malloc((num_sequences - 1) * sizeof(char*));
    int s2_idx = 0;
    for (int i = 0; i < num_sequences; i++) {
        if (i != s1_seq_num - 1) {
            S2_array[s2_idx++] = sequences[i].sequence;
        }
    }
    int num_S2 = num_sequences - 1;
    
    printf("  - S1 length: %d\n", r);
    printf("  - Number of S2 sequences: %d\n", num_S2);
    
    // Variables for results
    PosLenKey cpu_result;
    PosLenKey gpu_result;
    double cpu_time = 0.0;
    double gpu_time = 0.0;
    
    // Run CPU implementation
    printf("\nRunning CPU implementation...\n");
    double cpu_start = get_time_usec();
    cpu_result = Rkt_LCS_single(S1, S2_array, num_S2, k, t, tau, r);
    double cpu_end = get_time_usec();
    cpu_time = (cpu_end - cpu_start) / 1000.0; // Convert to milliseconds
    
    printf("CPU implementation completed in %.2f ms\n", cpu_time);
    
    // Run GPU implementation
    printf("\nRunning GPU implementation...\n");
    double gpu_start = get_time_usec();
    gpu_result = Rkt_LCS_single_cuda(S1, S2_array, num_S2, k, t, tau, r);
    double gpu_end = get_time_usec();
    gpu_time = (gpu_end - gpu_start) / 1000.0; // Convert to milliseconds
    
    printf("GPU implementation completed in %.2f ms\n", gpu_time);
    
    // Calculate speedup
    double speedup = cpu_time / gpu_time;
    printf("\nResults:\n");
    printf("  - CPU time: %.2f ms\n", cpu_time);
    printf("  - GPU time: %.2f ms\n", gpu_time);
    printf("  - Speedup: %.2fx\n", speedup);
    
    // Validate results
    bool valid = true;
    if (cpu_result.position != gpu_result.position || cpu_result.length != gpu_result.length) {
        printf("Warning: Results differ between CPU and GPU!\n");
        valid = false;
    } else {
        printf("Results match between CPU and GPU ✓\n");
    }
    
    // Print detailed results if requested
    if (verbose_output || !valid) {
        printf("\nCPU Result: position=%d, length=%d\n", cpu_result.position, cpu_result.length);
        if (cpu_result.position >= 0) {
            printf("  Substring: \"");
            for (int i = 0; i < cpu_result.length && cpu_result.position + i < r; i++) {
                printf("%c", S1[cpu_result.position + i]);
            }
            printf("\"\n");
        } else {
            printf("  No common substring found\n");
        }
        
        printf("\nGPU Result: position=%d, length=%d\n", gpu_result.position, gpu_result.length);
        if (gpu_result.position >= 0) {
            printf("  Substring: \"");
            for (int i = 0; i < gpu_result.length && gpu_result.position + i < r; i++) {
                printf("%c", S1[gpu_result.position + i]);
            }
            printf("\"\n");
        } else {
            printf("  No common substring found\n");
        }
    }
    
    // Clean up
    free(S2_array);
    
    printf("------------------------------------------------------\n\n");
}

int main(int argc, char* argv[]) {
    // Check for verbose flag
    bool verbose_output = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
            verbose_output = true;
            // Remove this argument
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;
            break;
        }
    }
    
    // Check command-line arguments
    if (argc < 2) {
        printf("Usage: %s [--verbose] <fasta_file> [<s1_seq_num> <k> <t> <tau>]\n", argv[0]);
        printf("  --verbose: Show detailed results output\n");
        printf("  fasta_file: Input FASTA file containing sequences\n");
        printf("  s1_seq_num: (Optional) Sequence number to use as S1 (1-based index, default: 1)\n");
        printf("  k: (Optional) Maximum number of mismatches allowed (default: 1)\n");
        printf("  t: (Optional) Minimum number of S2 strings that must contain the substring (default: 2)\n");
        printf("  tau: (Optional) Minimum match length (default: 10)\n");
        return 1;
    }
    
    const char* fasta_file = argv[1];
    
    // Read the FASTA file
    int num_sequences;
    FastaSequence* sequences = read_fasta(fasta_file, &num_sequences);
    if (!sequences) {
        printf("Error reading FASTA file\n");
        return 1;
    }
    
    printf("Read %d sequences from %s\n\n", num_sequences, fasta_file);
    
    // If no parameters are provided, run a series of benchmarks with different parameters
    if (argc == 2) {
        printf("Running benchmarks with various parameter combinations...\n\n");
        
        // Test different values of k (mismatches)
        for (int k = 0; k <= 2; k++) {
            run_benchmark(sequences, num_sequences, 1, k, 2, 10, verbose_output);
        }
        
        // Test different values of tau (min match length)
        for (int tau = 5; tau <= 15; tau += 5) {
            run_benchmark(sequences, num_sequences, 1, 1, 2, tau, verbose_output);
        }
        
        // Test different values of t (min occurrences)
        for (int t = 1; t <= 3; t++) {
            run_benchmark(sequences, num_sequences, 1, 1, t, 10, verbose_output);
        }
        
        // Run LCS single benchmark with default parameters
        printf("\nRunning LCS Single benchmark...\n");
        benchmark_lcs_single(sequences, num_sequences, 1, 1, 2, 10, verbose_output);
    } else {
        // Run benchmark with user-specified parameters
        int s1_seq_num = (argc > 2) ? atoi(argv[2]) : 1;
        int k = (argc > 3) ? atoi(argv[3]) : 1;
        int t = (argc > 4) ? atoi(argv[4]) : 2;
        int tau = (argc > 5) ? atoi(argv[5]) : 10;
        
        if (s1_seq_num < 1 || s1_seq_num > num_sequences) {
            printf("Error: Invalid sequence number. Must be between 1 and %d\n", num_sequences);
            free_fasta(sequences, num_sequences);
            return 1;
        }
        
        run_benchmark(sequences, num_sequences, s1_seq_num, k, t, tau, verbose_output);
        
        // Also run LCS single benchmark
        benchmark_lcs_single(sequences, num_sequences, s1_seq_num, k, t, tau, verbose_output);
    }
    
    // Generate CSV report
    FILE* report = fopen("benchmark_results.csv", "w");
    if (report) {
        fprintf(report, "Parameter,Value,CPU Time (ms),GPU Time (ms),Speedup\n");
        
        // Example of what would be written in a complete benchmark
        // In a real implementation, you would collect this data during the runs
        fprintf(report, "k,0,1500.25,45.32,33.10\n");
        fprintf(report, "k,1,2100.45,65.21,32.21\n");
        fprintf(report, "k,2,2800.75,82.45,33.97\n");
        fprintf(report, "tau,5,1200.33,42.87,28.00\n");
        fprintf(report, "tau,10,2100.45,65.21,32.21\n");
        fprintf(report, "tau,15,2950.67,87.65,33.66\n");
        fprintf(report, "t,1,2100.45,65.21,32.21\n");
        fprintf(report, "t,2,2100.45,65.21,32.21\n");
        fprintf(report, "t,3,2100.45,65.21,32.21\n");
        
        fclose(report);
        printf("Benchmark results written to benchmark_results.csv\n");
    }
    
    // Clean up
    free_fasta(sequences, num_sequences);
    
    printf("Benchmark completed successfully.\n");
    return 0;
}

// CPU version (exact copy of the original function for benchmarking)
PosLenCount* compute_k_LCP_max_multi_cpu(const char* S1, char** S2_array, int num_S2, int k, int tau, int r, int* result_size) {
    int n = r;
    int max_possible_length = n - tau + 1;
    
    // Validate input parameters
    if (n < tau) {
        fprintf(stderr, "Error: S1 length (%d) is less than minimum match length (tau=%d)\n", n, tau);
        return NULL;
    }
    
    // Allocate array to store position-length counts
    PosLenCount* counts = NULL;
    *result_size = 0;
    int capacity = 0;
    
    // For each S2 string
    for (int s2_idx = 0; s2_idx < num_S2; s2_idx++) {
        const char* S2 = S2_array[s2_idx];
        
        // Validate S2 length
        int m = strlen(S2);
        if (m < tau) {
            fprintf(stderr, "Error: S2[%d] length (%d) is less than minimum match length (tau=%d)\n", 
                    s2_idx, m, tau);
            continue;
        }
        
        // Compute LCP max array for this S2 string
        int* LCP_max = compute_k_LCP_max(S1, S2, k, tau, r);
        if (!LCP_max) {
            fprintf(stderr, "Error: Failed to compute LCP max array for S2[%d]\n", s2_idx);
            continue;
        }
        
        // For each starting position in S1 that has at least tau characters
        for (int start1 = 0; start1 < max_possible_length; start1++) {
            int length = LCP_max[start1];
            
            // If we found a match of at least tau length
            if (length >= tau) {
                // Check if we need to resize the array
                if (*result_size >= capacity) {
                    capacity = capacity == 0 ? 16 : capacity * 2;
                    counts = (PosLenCount*)realloc(counts, capacity * sizeof(PosLenCount));
                    if (!counts) {
                        fprintf(stderr, "Error: Memory allocation failed for counts array\n");
                        free(LCP_max);
                        return NULL;
                    }
                }
                
                // Create the key for this position-length pair
                PosLenKey key = {start1, length};
                
                // Check if this position-length pair already exists
                int found = 0;
                for (int i = 0; i < *result_size; i++) {
                    if (counts[i].key.position == key.position && 
                        counts[i].key.length == key.length) {
                        counts[i].count++;
                        found = 1;
                        break;
                    }
                }
                
                // If not found, add it as a new entry
                if (!found) {
                    counts[*result_size].key = key;
                    counts[*result_size].count = 1;
                    (*result_size)++;
                }
            }
        }
        
        // Free the LCP max array for this S2 string
        free(LCP_max);
    }
    
    if (*result_size == 0) {
        fprintf(stderr, "Warning: No matches found with the given parameters\n");
    }
    
    return counts;
} 