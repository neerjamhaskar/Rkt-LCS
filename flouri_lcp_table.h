#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Queue structure and operations remain the same
typedef struct {
    int* elements;
    int capacity;
    int size;
    int front;
    int rear;
} Queue;

Queue* createQueue(int capacity) {
    Queue* queue = (Queue*)malloc(sizeof(Queue));
    queue->elements = (int*)malloc(capacity * sizeof(int));
    queue->capacity = capacity;
    queue->size = 0;
    queue->front = 0;
    queue->rear = -1;
    return queue;
}

void enqueue(Queue* Q, int value) {
    if (Q->size == Q->capacity) return;
    Q->rear = (Q->rear + 1) % Q->capacity;
    Q->elements[Q->rear] = value;
    Q->size++;
}

void dequeue(Queue* Q) {
    if (Q->size == 0) return;
    Q->front = (Q->front + 1) % Q->capacity;
    Q->size--;
}

int minQueue(Queue* Q) {
    if (Q->size == 0) return -1;
    int min = Q->elements[Q->front];
    int count = Q->size;
    int index = Q->front;
    while (count > 0) {
        if (Q->elements[index] < min) {
            min = Q->elements[index];
        }
        index = (index + 1) % Q->capacity;
        count--;
    }
    return min;
}

// helper function 
void print2DArray(int** arr, int rows, int cols);

int** compute_k_LCP(const char* S1, const char* S2, int k) {
    int n = strlen(S1);
    int m = strlen(S2);

    // Allocate and initialize LCP table
    int** LCP = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        LCP[i] = (int*)malloc(m * sizeof(int));
    }

    Queue* Q = createQueue(k);

    // Compute LCP for each suffix pair
    for (int start1 = 0; start1 < n; start1++) {
        for (int start2 = 0; start2 < m; start2++) {
            // Clear queue
            while (Q->size > 0) {
                dequeue(Q);
            }

            int p = 0;
            int max_length = 0;

            // Compare characters and count mismatches
            while ((start1 + p < n) && (start2 + p < m)) {
                if (S1[start1 + p] != S2[start2 + p]) {
                    if (Q->size == k) {
                        break;  // Stop when we exceed k mismatches
                    }
                    enqueue(Q, p);
                }
                p++;
                max_length = p;
            }

            LCP[start1][start2] = max_length;
        }
    }

    // Clean up
    free(Q->elements);
    free(Q);

    return LCP;  // Return the LCP table
}

// LengthStats[L, j]
int** compute_LengthStat(int*** LCP_i, char** S, int m, int p, int i) {
    const int li = strlen(S[i]);
    const int L = li - p;     // Length of the longest common prefix. 1~L

    // Allocate and zero-initialize the array of pointers
    int** LengthStat = (int**)calloc(L, sizeof(int*));
    for (int l = 1; l <= L; l++)
        LengthStat[l - 1] = (int*)calloc(m + 1, sizeof(int));

    // Initialize
    for (int j = i; j < m; j++) {
        const char* sj = S[j];
        const int lj = strlen(sj);
        for (int q = 0; q < lj; q++) {
            const int l = LCP_i[j][p][q];
            if (l > 0)
                LengthStat[l - 1][j] = 1;
        }
    }

    // Compute
    int sum = 0;
    for (int j = 0; j < m; j++) {
        sum += LengthStat[L - 1][j];
    }
    LengthStat[L - 1][m] = sum;

    for (int l = L - 1; l > 0; l--) {
        int sum = 0;
        for (int j = i; j < m; j++) {
            LengthStat[l - 1][j] |= LengthStat[l][j];
            sum += LengthStat[l - 1][j];
        }
        LengthStat[l - 1][m] = sum;
    }

    return LengthStat;  // Return the LCP table
}

/* 
 * Returns a new string containing a substring of 'str' starting at offset 'p' with length 'l'.
 * If 'p' or 'l' are out of bounds, the behavior is undefined (you might add error-checking as needed).
 */
char * substring(const char *str, int p, int l) {
    // Allocate memory for l characters plus the null terminator.
    char *result = (char *)malloc(l + 1);
    if (result == NULL) {
        // Memory allocation failed.
        return NULL;
    }
    
    // Copy l characters starting at (str + p) into result.
    strncpy(result, str + p, l);
    
    // Ensure the string is null-terminated.
    result[l] = '\0';
    
    return result;
}

char * Rkt_LCS(char* S[], int m, int k, int t) {
    int len_max = 0;
    char * result = NULL;
    for (int i = 0; i < m; i++) {
        // anchor on S[i]
        // compute LCP_i table
        int*** LCP_i = (int***) malloc(m * sizeof(int**));  // but only i...m-1 has table
        for (int j = i; j < m; j++) {
            LCP_i[j] = compute_k_LCP(S[i], S[j], k);
        }

        int li = strlen(S[i]);
        for (int p = 0; p < li; p++) {
            // compute_LengthStat
            int** LengthStat = compute_LengthStat(LCP_i, S, m, p, i);

            // check if LengthStat[L, m] >= t
            for (int l = li - p; l >= 1; l--) {  
                if (LengthStat[l - 1][m] >= t && l > len_max) {
                    // printf("Qualified string found at %d: %.*s with length %d\n", p, l, S[i] + p, l);
                    len_max = l;
                    result = substring(S[i], p, l);
                }
            }
        }
    }
    return result;
}



