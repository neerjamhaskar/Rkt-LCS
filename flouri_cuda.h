#ifndef FLOURI_CUDA_H
#define FLOURI_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

// Function declarations
int** compute_k_lcp_cuda(const char* S1, const char* S2, int k);
void free_2d_array(int** arr, int rows);

#ifdef __cplusplus
}
#endif

#endif // FLOURI_CUDA_H 