#ifndef FLOURI_CPU_H
#define FLOURI_CPU_H

#ifdef __cplusplus
extern "C" {
#endif

// Function declarations
int** compute_k_LCP(const char* S1, const char* S2, int k);
void free_2d_array(int** arr, int rows);

#ifdef __cplusplus
}
#endif

#endif // FLOURI_CPU_H 