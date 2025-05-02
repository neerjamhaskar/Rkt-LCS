#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "flouri_cuda.h"

// Error checking macro for CUDA calls
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(error)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

// Queue structure for GPU with fixed-size array
struct Queue {
    int elements[32];  // Fixed size array, assuming k <= 32
    int capacity;
    int size;
    int front;
    int rear;
};

// CUDA kernel for computing k-difference LCP table
__global__ void compute_k_lcp_kernel(
    const char* S1,
    const char* S2,
    int* LCP,
    int n,
    int m,
    int k
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i >= n || j >= m) return;

    // Create queue for this thread
    Queue Q;
    Q.capacity = k;
    Q.size = 0;
    Q.front = 0;
    Q.rear = -1;

    int p = 0;
    int max_length = 0;

    // Compare characters and count mismatches
    while ((i + p < n) && (j + p < m)) {
        if (S1[i + p] != S2[j + p]) {
            if (Q.size == k) {
                break;  // Stop when we exceed k mismatches
            }
            // Enqueue
            Q.rear = (Q.rear + 1) % k;
            Q.elements[Q.rear] = p;
            Q.size++;
        }
        p++;
        max_length = p;
    }

    // Store result in flattened array
    LCP[i * m + j] = max_length;
}

// Main function to compute k-difference LCP table using CUDA
extern "C" int** compute_k_lcp_cuda(const char* S1, const char* S2, int k) {
    int n = strlen(S1);
    int m = strlen(S2);

    // Allocate device memory
    char* d_S1;
    char* d_S2;
    int* d_LCP;
    int** h_LCP;

    // Use pinned memory for faster transfers
    char* h_S1_pinned, *h_S2_pinned;

    // Calculate total memory needed
    size_t total_memory_needed = (n * sizeof(char)) + (m * sizeof(char)) + (n * m * sizeof(int));
    size_t free_memory, total_memory;
    cudaMemGetInfo(&free_memory, &total_memory);
    
    if (total_memory_needed > free_memory) {
        fprintf(stderr, "Not enough GPU memory. Need %zu bytes, but only %zu available\n", 
                total_memory_needed, free_memory);
        return NULL;
    }

    CUDA_CHECK(cudaMallocHost(&h_S1_pinned, n * sizeof(char)));
    CUDA_CHECK(cudaMallocHost(&h_S2_pinned, m * sizeof(char)));
    memcpy(h_S1_pinned, S1, n * sizeof(char));
    memcpy(h_S2_pinned, S2, m * sizeof(char));

    // Allocate device memory with error checking
    cudaError_t err;
    err = cudaMalloc(&d_S1, n * sizeof(char));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate d_S1: %s\n", cudaGetErrorString(err));
        return NULL;
    }
    
    err = cudaMalloc(&d_S2, m * sizeof(char));
    if (err != cudaSuccess) {
        cudaFree(d_S1);
        fprintf(stderr, "Failed to allocate d_S2: %s\n", cudaGetErrorString(err));
        return NULL;
    }
    
    err = cudaMalloc(&d_LCP, n * m * sizeof(int));
    if (err != cudaSuccess) {
        cudaFree(d_S1);
        cudaFree(d_S2);
        fprintf(stderr, "Failed to allocate d_LCP: %s\n", cudaGetErrorString(err));
        return NULL;
    }

    // Allocate host memory for result
    h_LCP = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        h_LCP[i] = (int*)malloc(m * sizeof(int));
    }

    // Copy input data to device using pinned memory
    CUDA_CHECK(cudaMemcpy(d_S1, h_S1_pinned, n * sizeof(char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_S2, h_S2_pinned, m * sizeof(char), cudaMemcpyHostToDevice));

    // Calculate grid and block dimensions
    dim3 block_size(16, 16);  // Smaller block size for better occupancy
    dim3 grid_size((n + block_size.x - 1) / block_size.x, 
                   (m + block_size.y - 1) / block_size.y);

    // Launch kernel
    compute_k_lcp_kernel<<<grid_size, block_size>>>(d_S1, d_S2, d_LCP, n, m, k);

    // Copy result back to host row by row
    for (int i = 0; i < n; i++) {
        CUDA_CHECK(cudaMemcpy(h_LCP[i], d_LCP + i * m, m * sizeof(int), cudaMemcpyDeviceToHost));
    }

    // Free device memory
    CUDA_CHECK(cudaFree(d_S1));
    CUDA_CHECK(cudaFree(d_S2));
    CUDA_CHECK(cudaFree(d_LCP));
    CUDA_CHECK(cudaFreeHost(h_S1_pinned));
    CUDA_CHECK(cudaFreeHost(h_S2_pinned));

    return h_LCP;
}

// Function to free 2D array
extern "C" void free_2d_array(int** arr, int rows) {
    for (int i = 0; i < rows; i++) {
        free(arr[i]);
    }
    free(arr);
}

// int main(int argc, char* argv[]) {
//     if (argc != 3) {
//         printf("Usage: %s <string1> <string2>\n", argv[0]);
//         return 1;
//     }

//     const char* str1 = argv[1];
//     const char* str2 = argv[2];
//     int k = 2;  // Default k value
    
//     printf("String 1: %s\n", str1);
//     printf("String 2: %s\n", str2);
//     printf("k: %d\n", k);
    
//     // Compute k-difference LCP table
//     int** lcp_table = compute_k_lcp_cuda(str1, str2, k);
    
//     // Print the LCP table
//     printf("\nLCP Table (k=%d):\n", k);
//     printf("   ");
//     for (int j = 0; j < strlen(str2); j++) {
//         printf("%3d ", j);
//     }
//     printf("\n");

//     for (int i = 0; i < strlen(str1); i++) {
//         printf("%2d: ", i);
//         for (int j = 0; j < strlen(str2); j++) {
//             printf("%3d ", lcp_table[i][j]);
//         }
//         printf("\n");
//     }
    
//     // Free memory
//     free_2d_array(lcp_table, strlen(str1));
    
//     return 0;
// } 