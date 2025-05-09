#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include "flouri_cpu.h"

// Constants for GPU processing
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error in %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

// CUDA kernel for computing LCP max with k mismatches
__global__ void compute_k_LCP_max_kernel(
    const char* S1,                  // S1 sequence
    const char* S2_buffer,           // Flattened buffer of all S2 sequences
    const int* S2_offsets,           // Starting offset of each S2 in the buffer
    const int* S2_lengths,           // Length of each S2 sequence
    int num_S2,                      // Number of S2 sequences
    int S1_length,                   // Length of S1
    int k,                           // Max mismatches allowed
    int tau,                         // Minimum match length
    int* results,                    // Output buffer: (S2_idx, start1, length)
    int* result_count                // Atomic counter for results
) {
    // Calculate the global thread index
    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Each thread processes one (S2_idx, start1) pair
    int S2_idx = global_idx / (S1_length - tau + 1);
    int start1 = global_idx % (S1_length - tau + 1);
    
    // Check if this is a valid thread
    if (S2_idx >= num_S2)
        return;
    
    // Get details for this S2 sequence
    int S2_offset = S2_offsets[S2_idx];
    int S2_length = S2_lengths[S2_idx];
    
    // Skip if S2 is too short
    if (S2_length < tau)
        return;
    
    // Keep track of mismatches
    int mismatches[64];  // Preallocated array for mismatch positions (limited to 64)
    int mismatch_count = 0;
    
    // For each possible starting position in S2, find longest match
    int longest = 0;
    
    for (int start2 = 0; start2 < S2_length; start2++) {
        // Reset mismatch count for this start2 position
        mismatch_count = 0;
        
        // Phase 1: Process the first tau characters
        int p;
        for (p = 0; p < tau; p++) {
            if (start1 + p >= S1_length || start2 + p >= S2_length)
                break;  // End of one of the strings
                
            if (S1[start1 + p] != S2_buffer[S2_offset + start2 + p]) {
                if (mismatch_count == k)
                    break;  // Exceeded allowed mismatches
                    
                if (mismatch_count < 64)  // Guard against buffer overflow
                    mismatches[mismatch_count++] = p;
            }
        }
        
        // If we couldn't process at least tau characters, discard this pair
        if (p < tau) {
            continue;
        }
        
        // Phase 2: Continue comparing beyond the first tau characters
        while ((start1 + p < S1_length) && (start2 + p < S2_length)) {
            if (S1[start1 + p] != S2_buffer[S2_offset + start2 + p]) {
                if (mismatch_count == k)
                    break;  // Exceeded allowed mismatches
                    
                if (mismatch_count < 64)  // Guard against buffer overflow
                    mismatches[mismatch_count++] = p;
            }
            p++;
        }
        
        // Update longest match if current match is longer
        if (p > longest) {
            longest = p;
        }
    }
    
    // If we found a match of at least tau length
    if (longest >= tau) {
        // Atomic add to get the next available index in the results array
        int idx = atomicAdd(result_count, 1);
        
        // Store the result: S2_idx, start position, and length
        results[idx * 3] = S2_idx;
        results[idx * 3 + 1] = start1;
        results[idx * 3 + 2] = longest;
    }
}

// Host function to compute MaxLCP_k multi on GPU
// We'll export this function with extern "C" to avoid name mangling
extern "C" PosLenCount* compute_k_LCP_max_multi_cuda(
    const char* S1, 
    char** S2_array, 
    int num_S2, 
    int k, 
    int tau, 
    int r, 
    int* result_size
) {
    int S1_length = r;
    int max_possible_length = S1_length - tau + 1;
    
    // Validate input parameters
    if (S1_length < tau) {
        fprintf(stderr, "Error: S1 length (%d) is less than minimum match length (tau=%d)\n", 
                S1_length, tau);
        return NULL;
    }
    
    // Report initial progress
    printf("Computing LCP results (0%%)\r");
    fflush(stdout);
    
    // Calculate total size needed for S2 buffer and prepare offsets and lengths
    int* S2_lengths = (int*)malloc(num_S2 * sizeof(int));
    int* S2_offsets = (int*)malloc(num_S2 * sizeof(int));
    int total_S2_length = 0;
    
    for (int i = 0; i < num_S2; i++) {
        S2_lengths[i] = strlen(S2_array[i]);
        S2_offsets[i] = total_S2_length;
        total_S2_length += S2_lengths[i];
    }
    
    // Flatten S2 sequences into a single buffer
    char* S2_buffer = (char*)malloc(total_S2_length * sizeof(char));
    for (int i = 0; i < num_S2; i++) {
        memcpy(S2_buffer + S2_offsets[i], S2_array[i], S2_lengths[i]);
    }
    
    // Allocate device memory
    char* d_S1;
    char* d_S2_buffer;
    int* d_S2_offsets;
    int* d_S2_lengths;
    int* d_results;  // Format: [S2_idx, start1, length, S2_idx, start1, length, ...]
    int* d_result_count;
    
    // Estimate max possible results (worst case: every position in every S2 matches)
    int max_results = num_S2 * max_possible_length;
    
    // Allocate device memory
    CUDA_CHECK(cudaMalloc((void**)&d_S1, S1_length * sizeof(char)));
    CUDA_CHECK(cudaMalloc((void**)&d_S2_buffer, total_S2_length * sizeof(char)));
    CUDA_CHECK(cudaMalloc((void**)&d_S2_offsets, num_S2 * sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_S2_lengths, num_S2 * sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, max_results * 3 * sizeof(int)));  // 3 integers per result
    CUDA_CHECK(cudaMalloc((void**)&d_result_count, sizeof(int)));
    
    // Initialize result count to 0
    int h_result_count = 0;
    CUDA_CHECK(cudaMemcpy(d_result_count, &h_result_count, sizeof(int), cudaMemcpyHostToDevice));
    
    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_S1, S1, S1_length * sizeof(char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_S2_buffer, S2_buffer, total_S2_length * sizeof(char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_S2_offsets, S2_offsets, num_S2 * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_S2_lengths, S2_lengths, num_S2 * sizeof(int), cudaMemcpyHostToDevice));
    
    // Calculate kernel launch parameters
    int total_threads = num_S2 * max_possible_length;
    int threads_per_block = 256;
    int blocks = (total_threads + threads_per_block - 1) / threads_per_block;
    
    printf("Launching CUDA kernel with %d blocks, %d threads per block (%d total threads)\n", 
           blocks, threads_per_block, total_threads);
    
    // Launch kernel
    compute_k_LCP_max_kernel<<<blocks, threads_per_block>>>(
        d_S1, d_S2_buffer, d_S2_offsets, d_S2_lengths,
        num_S2, S1_length, k, tau, d_results, d_result_count
    );
    
    // Check for kernel launch errors
    CUDA_CHECK(cudaGetLastError());
    
    // Wait for kernel to finish
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Copy result count back to host
    CUDA_CHECK(cudaMemcpy(&h_result_count, d_result_count, sizeof(int), cudaMemcpyDeviceToHost));
    
    // Update progress to show completion
    printf("Computing LCP results (100%%)\n");
    
    printf("CUDA kernel produced %d raw results\n", h_result_count);
    
    // Check if we have valid results
    if (h_result_count == 0) {
        fprintf(stderr, "Warning: No matches found with the given parameters\n");
        
        // Clean up
        cudaFree(d_S1);
        cudaFree(d_S2_buffer);
        cudaFree(d_S2_offsets);
        cudaFree(d_S2_lengths);
        cudaFree(d_results);
        cudaFree(d_result_count);
        
        free(S2_buffer);
        free(S2_offsets);
        free(S2_lengths);
        
        *result_size = 0;
        return NULL;
    }
    
    // Allocate memory for results on host
    int* h_results = (int*)malloc(h_result_count * 3 * sizeof(int));
    
    // Copy results back to host
    CUDA_CHECK(cudaMemcpy(h_results, d_results, h_result_count * 3 * sizeof(int), cudaMemcpyDeviceToHost));
    
    // Free device memory
    cudaFree(d_S1);
    cudaFree(d_S2_buffer);
    cudaFree(d_S2_offsets);
    cudaFree(d_S2_lengths);
    cudaFree(d_results);
    cudaFree(d_result_count);
    
    // Free temporary host buffers
    free(S2_buffer);
    free(S2_offsets);
    free(S2_lengths);
    
    // Process raw results into PosLenCount format (same as CPU version)
    PosLenCount* counts = NULL;
    *result_size = 0;
    int capacity = 0;
    
    for (int i = 0; i < h_result_count; i++) {
        int S2_idx = h_results[i * 3];
        int start1 = h_results[i * 3 + 1];
        int length = h_results[i * 3 + 2];
        // For all lengths from tau up to the maximum found
        for (int l = tau; l <= length; l++) {
            PosLenKey key = {start1, l};
            // Check if we need to resize the array
            if (*result_size >= capacity) {
                capacity = capacity == 0 ? 16 : capacity * 2;
                counts = (PosLenCount*)realloc(counts, capacity * sizeof(PosLenCount));
                if (!counts) {
                    fprintf(stderr, "Error: Memory allocation failed for counts array\n");
                    free(h_results);
                    return NULL;
                }
            }
            // Check if this position-length pair already exists
            int found = 0;
            for (int j = 0; j < *result_size; j++) {
                if (counts[j].key.position == key.position && 
                    counts[j].key.length == key.length) {
                    counts[j].count++;
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
    
    // Clean up
    free(h_results);
    
    return counts;
}

// CUDA-accelerated version of Rkt_LCS_single
//
// For a single string S1 and an array of S2 strings, this function finds
// the Rkt-LCS in S1 (candidate) that appears in at least t S2 strings
// (allowing up to k mismatches). It returns a PosLenKey containing the
// position and length of the longest such substring.
// This version leverages the GPU implementation for better performance.
extern "C" PosLenKey Rkt_LCS_single_cuda(const char* S1, char** S2_array, int num_S2, int k, int t, int tau, int r) {
    int result_size;
    PosLenCount* results = compute_k_LCP_max_multi_cuda(S1, S2_array, num_S2, k, tau, r, &result_size);
    
    if (!results) {
        fprintf(stderr, "Failed to compute multi-string LCP results on GPU\n");
        return (PosLenKey){-1, -1};  // Return invalid key to indicate error
    }
    
    // Find the longest substring that appears in at least t strings
    PosLenKey best_key = {-1, -1};  // Initialize with invalid values
    int max_length = 0;
    
    for (int i = 0; i < result_size; i++) {
        if (results[i].count >= t && results[i].key.length > max_length) {
            max_length = results[i].key.length;
            best_key = results[i].key;
        }
    }
    
    // Clean up
    free(results);
    
    return best_key;
}
