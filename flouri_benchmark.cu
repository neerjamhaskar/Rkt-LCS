#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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
int** compute_k_lcp_cuda(const char* S1, const char* S2, int k, double* transfer_time, double* kernel_time) {
    int n = strlen(S1);
    int m = strlen(S2);

    // Allocate device memory
    char* d_S1;
    char* d_S2;
    int* d_LCP;
    int** h_LCP;

    // Start timing memory transfers
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start));

    // Use pinned memory for faster transfers
    char* h_S1_pinned, *h_S2_pinned;
    CUDA_CHECK(cudaMallocHost(&h_S1_pinned, n * sizeof(char)));
    CUDA_CHECK(cudaMallocHost(&h_S2_pinned, m * sizeof(char)));
    memcpy(h_S1_pinned, S1, n * sizeof(char));
    memcpy(h_S2_pinned, S2, m * sizeof(char));

    CUDA_CHECK(cudaMalloc(&d_S1, n * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&d_S2, m * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&d_LCP, n * m * sizeof(int)));

    // Allocate host memory for result
    h_LCP = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        h_LCP[i] = (int*)malloc(m * sizeof(int));
    }

    // Copy input data to device using pinned memory
    CUDA_CHECK(cudaMemcpy(d_S1, h_S1_pinned, n * sizeof(char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_S2, h_S2_pinned, m * sizeof(char), cudaMemcpyHostToDevice));

    // Record transfer time
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float transfer_ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&transfer_ms, start, stop));
    *transfer_time = transfer_ms / 1000.0;

    // Calculate grid and block dimensions
    dim3 block_size(16, 16);  // Smaller block size for better occupancy
    dim3 grid_size((n + block_size.x - 1) / block_size.x, 
                   (m + block_size.y - 1) / block_size.y);

    // Start timing kernel execution
    CUDA_CHECK(cudaEventRecord(start));

    // Launch kernel
    compute_k_lcp_kernel<<<grid_size, block_size>>>(d_S1, d_S2, d_LCP, n, m, k);

    // Record kernel time
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float kernel_ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&kernel_ms, start, stop));
    *kernel_time = kernel_ms / 1000.0;

    // Start timing result transfer
    CUDA_CHECK(cudaEventRecord(start));

    // Copy result back to host row by row
    for (int i = 0; i < n; i++) {
        CUDA_CHECK(cudaMemcpy(h_LCP[i], d_LCP + i * m, m * sizeof(int), cudaMemcpyDeviceToHost));
    }

    // Record result transfer time
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float result_transfer_ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&result_transfer_ms, start, stop));
    *transfer_time += result_transfer_ms / 1000.0;

    // Free device memory
    CUDA_CHECK(cudaFree(d_S1));
    CUDA_CHECK(cudaFree(d_S2));
    CUDA_CHECK(cudaFree(d_LCP));
    CUDA_CHECK(cudaFreeHost(h_S1_pinned));
    CUDA_CHECK(cudaFreeHost(h_S2_pinned));

    // Clean up CUDA events
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    return h_LCP;
}

// CPU implementation of k-difference LCP table
int** compute_k_lcp_cpu(const char* S1, const char* S2, int k) {
    int n = strlen(S1);
    int m = strlen(S2);
    
    // Allocate result array
    int** LCP = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        LCP[i] = (int*)malloc(m * sizeof(int));
    }

    // Compute LCP table
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            Queue Q;
            Q.capacity = k;
            Q.size = 0;
            Q.front = 0;
            Q.rear = -1;

            int p = 0;
            int max_length = 0;

            while ((i + p < n) && (j + p < m)) {
                if (S1[i + p] != S2[j + p]) {
                    if (Q.size == k) {
                        break;
                    }
                    Q.rear = (Q.rear + 1) % k;
                    Q.elements[Q.rear] = p;
                    Q.size++;
                }
                p++;
                max_length = p;
            }

            LCP[i][j] = max_length;
        }
    }

    return LCP;
}

// Function to generate random string
char* generate_random_string(int length) {
    char* str = (char*)malloc((length + 1) * sizeof(char));
    const char charset[] = "abcdefghijklmnopqrstuvwxyz";
    for (int i = 0; i < length; i++) {
        str[i] = charset[rand() % (sizeof(charset) - 1)];
    }
    str[length] = '\0';
    return str;
}

// Function to free 2D array
void free_2d_array(int** arr, int rows) {
    for (int i = 0; i < rows; i++) {
        free(arr[i]);
    }
    free(arr);
}

int main(int argc, char* argv[]) {
    // Initialize random seed
    srand(time(NULL));

    // Generate two random strings of length 5000 (increased for better GPU utilization)
    int str_length = 5000;
    char* str1 = generate_random_string(str_length);
    char* str2 = generate_random_string(str_length);
    int k = 2;

    printf("Benchmarking with strings of length %d\n", str_length);
    printf("k = %d\n", k);

    // CPU Implementation
    clock_t cpu_start = clock();
    int** cpu_result = compute_k_lcp_cpu(str1, str2, k);
    clock_t cpu_end = clock();
    double cpu_time = ((double)(cpu_end - cpu_start)) / CLOCKS_PER_SEC;

    // CUDA Implementation with timing breakdown
    double transfer_time = 0, kernel_time = 0;
    int** cuda_result = compute_k_lcp_cuda(str1, str2, k, &transfer_time, &kernel_time);
    double total_cuda_time = transfer_time + kernel_time;

    // Verify results match
    bool results_match = true;
    for (int i = 0; i < str_length; i++) {
        for (int j = 0; j < str_length; j++) {
            if (cpu_result[i][j] != cuda_result[i][j]) {
                results_match = false;
                break;
            }
        }
        if (!results_match) break;
    }

    // Print results
    printf("\nPerformance Results:\n");
    printf("CPU Time: %.4f seconds\n", cpu_time);
    printf("CUDA Time Breakdown:\n");
    printf("  - Memory Transfer: %.4f seconds\n", transfer_time);
    printf("  - Kernel Execution: %.4f seconds\n", kernel_time);
    printf("  - Total CUDA Time: %.4f seconds\n", total_cuda_time);
    printf("Speedup: %.2fx\n", cpu_time / total_cuda_time);
    printf("Results Match: %s\n", results_match ? "Yes" : "No");

    // Free memory
    free(str1);
    free(str2);
    free_2d_array(cpu_result, str_length);
    free_2d_array(cuda_result, str_length);

    return 0;
} 