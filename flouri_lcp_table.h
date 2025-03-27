#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//#region <Helper DS and functions>
// Helper DS: Queue
//
// Queue structure and operations
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

// Helper function: print2DArray
//
// Prints a 2D integer array with the specified number of rows and columns.
void print2DArray(int** arr, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        printf("Row %d: ", i);
        for (int j = 0; j < cols; j++) {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }
}

// Helper function: printArray
//
// Prints an integer array with the specified number of items.
void printArray(int* arr, int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

// Helper function: substring
//
// Returns a newly allocated substring of 'str' starting at offset p with length l.
char * substring(const char *str, int p, int l) {
    char *result = (char *)malloc(l + 1);
    if (!result) return NULL;
    strncpy(result, str + p, l);
    result[l] = '\0';
    return result;
}

// Helper function: free2DArray
//
// Frees a 2D array that was allocated using dynamically.
void free2DArray(int** arr, int rows) {
    for (int i = 0; i < rows; i++) {
        free(arr[i]);
    }
    free(arr);
}

//#endregion




/***  Main functions  ***/
/*---------------------------------------------------------*/

// compute_k_LCP
//
// Computes a longest-common-prefix (LCP) table between S1 and S2,
// allowing up to k mismatches. The LCP table is of size
// strlen(S1) x strlen(S2). A temporary queue is used to record mismatch positions.
int** compute_k_LCP(const char* S1, const char* S2, int k) {
    int n = strlen(S1);
    int m = strlen(S2);

    // Allocate LCP table: n rows, each with m integers.
    int** LCP = (int**)malloc(n * sizeof(int*));
    if (!LCP) {
        fprintf(stderr, "Memory allocation failed for LCP table\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++) {
        LCP[i] = (int*)malloc(m * sizeof(int));
        if (!LCP[i]) {
            fprintf(stderr, "Memory allocation failed for LCP[%d]\n", i);
            for (int j = 0; j < i; j++) free(LCP[j]);
            free(LCP);
            exit(EXIT_FAILURE);
        }
    }

    // Create a queue to track mismatch positions (capacity = k)
    Queue* Q = createQueue(k);

    // For each pair of starting positions in S1 and S2:
    for (int start1 = 0; start1 < n; start1++) {
        for (int start2 = 0; start2 < m; start2++) {
            /* Reset the queue quickly instead of dequeueing each element */
            Q->size = 0;
            Q->front = 0;
            Q->rear = -1;

            int p = 0, max_length = 0;
            while ((start1 + p < n) && (start2 + p < m)) {
                if (S1[start1 + p] != S2[start2 + p]) {
                    if (Q->size == k) {
                        break;  // Exceeded allowed mismatches.
                    }
                    enqueue(Q, p);
                }
                p++;
                max_length = p;
            }
            LCP[start1][start2] = max_length;
        }
    }

    free(Q->elements);
    free(Q);
    return LCP;
}

// compute_k_LCP_max (memory improved, with tau)
// Computes a longest-common-prefix (LCP) array between S1 and S2,
// allowing up to k mismatches. An extra parameter tau is used so that
// starting positions in S1 with less than tau characters remaining are discarded.
// Moreover, the first tau characters are compared normally (mismatches counted)
// and only if they are processed (i.e. p reaches tau) do we continue.
// Thus, the returned LCP length will include those tau characters.
int* compute_k_LCP_max(const char* S1, const char* S2, int k, int tau, int r) {
    int n = r;
    int m = r;

    // Only consider S1 positions that have at least tau characters remaining.
    int resultSize = n - tau + 1;
    int* LCP = (int*)malloc(resultSize * sizeof(int));
    if (!LCP) {
        fprintf(stderr, "Memory allocation failed for LCP table\n");
        exit(EXIT_FAILURE);
    }

    // Create a queue to track mismatch positions (capacity = k)
    Queue* Q = createQueue(k);

    // For each starting position in S1 that has at least tau characters.
    for (int start1 = 0; start1 < resultSize; start1++) {
        int longest = 0;
        for (int start2 = 0; start2 < m; start2++) {
            // Reset the queue for this start2.
            Q->size = 0;
            Q->front = 0;
            Q->rear = -1;

            int p = 0;

            // Phase 1: Process the first tau characters.
            // We compare the first tau characters and record any mismatches.
            for (p = 0; p < tau; p++) {
                if (start1 + p >= n || start2 + p >= m)
                    break;  // End of one of the strings.
                if (S1[start1 + p] != S2[start2 + p]) {
                    if (Q->size == k) {
                        break;  // Exceeded allowed mismatches.
                    }
                    enqueue(Q, p);
                }
            }
            // If we couldn't process at least tau characters, discard this pair.
            if (p < tau) {
                continue;
            }

            // Phase 2: Continue comparing beyond the first tau characters.
            while ((start1 + p < n) && (start2 + p < m)) {
                if (S1[start1 + p] != S2[start2 + p]) {
                    if (Q->size == k) {
                        break;  // Exceeded allowed mismatches.
                    }
                    enqueue(Q, p);
                }
                p++;
            }
            if (p > longest) {
                longest = p;
            }
        }
        LCP[start1] = longest;
    }

    free(Q->elements);
    free(Q);
    return LCP;
}

// compute_LengthStat
//
// Given a 2D LCP_i array (for strings j from i to m-1) and the array S,
// this function computes a LengthStat table for S[i] starting at offset p.
// The table has L rows (where L = strlen(S[i]) - p) and (m+1) columns;
// the last column is used to store cumulative sums.
// p needs to be <= len(S[i]) - tau
int** compute_LengthStat(int** LCP_i, char** S, int m, int p, int i, int r) {
    const int li = r;
    const int L = li - p;     // Maximum possible length for a substring from S[i] starting at p.

    // Allocate a table of L rows and (m+1) columns; initialize to zero.
    int** LengthStat = (int**)calloc(L, sizeof(int*));
    if (!LengthStat) {
        fprintf(stderr, "Memory allocation failed for LengthStat\n");
        exit(EXIT_FAILURE);
    }
    for (int l = 1; l <= L; l++) {
        LengthStat[l - 1] = (int*)calloc(m + 1, sizeof(int));
        if (!LengthStat[l - 1]) {
            fprintf(stderr, "Memory allocation failed for LengthStat row %d\n", l - 1);
            for (int j = 0; j < l - 1; j++) free(LengthStat[j]);
            free(LengthStat);
            exit(EXIT_FAILURE);
        }
    }

    // Initialization: for j from i to m-1, mark positions where LCP_i indicates a match.
    for (int j = i; j < m; j++) {
        const int l = LCP_i[j - i][p];    // LCP_i[j-i] has length len(S[i]) - tau + 1, so p needs to be <= len(S[i]) - tau.
        if (l > 0)
            LengthStat[l - 1][j] = 1;
    }

    // Compute cumulative sums in the last column.
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

// Rkt_LCS 
//
// For an array S of m strings, this function anchors on each S[i] and computes
// a longest common substring (LCS) among the strings (allowing up to k mismatches)
// that occurs in at least t strings. It returns a newly allocated substring.
char * Rkt_LCS(char* S[], int m, int k, int t, int tau, int r) {
    int len_max = 0;
    char * result = NULL;

    // for each string S[i]
    for (int i = 0; i < m; i++) {
        int li = strlen(S[i]);
        // Allocate LCP_i for strings j from i to m-1.
        const int numTables = m - i;
        int** LCP_i = (int**)malloc(numTables * sizeof(int*));
        if (!LCP_i) {
            fprintf(stderr, "Memory allocation failed for LCP_i\n");
            exit(EXIT_FAILURE);
        }
        for (int j = i; j < m; j++) {
            LCP_i[j - i] = compute_k_LCP_max(S[i], S[j], k, tau, r);
        }

        // For each possible starting position p in S[i]:
        for (int p = 0; p < li; p++) {
            int** LengthStat = compute_LengthStat(LCP_i, S, m, p, i, r);
            int L = li - p;  // Maximum possible substring length from S[i] starting at p.

            // Check for each candidate length (from longest down to 1)
            for (int l = L; l >= 1; l--) {  
                if (LengthStat[l - 1][m] >= t && l > len_max) {
                    if (result != NULL)
                        free(result);
                    len_max = l;
                    result = substring(S[i], p, l);
                }
            }

            // Free the LengthStat table.
            free2DArray(LengthStat, L);
        }

        // Free each LCP table.
        free2DArray(LCP_i, numTables);
    }

    return result;
}

/*** End of main functions ***/
/*---------------------------------------------------------*/



