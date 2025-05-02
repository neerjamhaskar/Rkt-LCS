#ifndef FLOURI_CPU_H
#define FLOURI_CPU_H

#ifdef __cplusplus
extern "C" {
#endif

// Queue structure
typedef struct {
    int* elements;
    int capacity;
    int size;
    int front;
    int rear;
} Queue;

// Position-length key structure
typedef struct {
    int position;
    int length;
} PosLenKey;

// Position-length count structure
typedef struct {
    PosLenKey key;
    int count;
} PosLenCount;

// FASTA sequence structure
typedef struct {
    char* name;
    char* sequence;
} FastaSequence;

// Queue operations
Queue* createQueue(int capacity);
void enqueue(Queue* Q, int value);
void dequeue(Queue* Q);
int minQueue(Queue* Q);

// Helper functions
void print2DArray(int** arr, int rows, int cols);
void printArray(int* arr, int n);
char* substring(const char* str, int p, int l);
void free2DArray(int** arr, int rows);

// Main functions
int** compute_k_LCP(const char* S1, const char* S2, int k);
int* compute_k_LCP_max(const char* S1, const char* S2, int k, int tau, int r);
void free_2d_array(int** arr, int rows);
PosLenCount* compute_k_LCP_max_multi(const char* S1, char** S2_array, int num_S2, int k, int tau, int r, int* result_size);
PosLenKey Rkt_LCS_single(const char* S1, char** S2_array, int num_S2, int k, int t, int tau, int r);
int flouri_main(int argc, char* argv[]);

// FASTA operations
FastaSequence* read_fasta(const char* filename, int* num_sequences);
void free_fasta(FastaSequence* sequences, int num_sequences);

// CUDA implementation (declared separately)
PosLenCount* compute_k_LCP_max_multi_cuda(const char* S1, char** S2_array, int num_S2, int k, int tau, int r, int* result_size);
PosLenKey Rkt_LCS_single_cuda(const char* S1, char** S2_array, int num_S2, int k, int t, int tau, int r);

#ifdef __cplusplus
}
#endif

#endif // FLOURI_CPU_H 