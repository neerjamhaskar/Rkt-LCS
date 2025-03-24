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

// Modified to compute k-difference LCP for all suffix pairs
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


// Test main function
// int main() {
//     // Test sequences
//     const char* S1 = "ssstring1";
//     const char* S2 = "string2";
//     int k = 2;  // Allow 2 mismatches

//     printf("Testing compute_k_LCP with:\n");
//     printf("S1: %s\n", S1);
//     printf("S2: %s\n", S2);
//     printf("k: %d\n\n", k);

//     // Compute LCP table
//     int** LCP = compute_k_LCP(S1, S2, k);
//     if (!LCP) {
//         printf("Failed to compute LCP table\n");
//         return 1;
//     }

//     // Print LCP table
//     printf("LCP Table (k=%d):\n", k);
//     printf("   ");
//     for (int j = 0; j < strlen(S2); j++) {
//         printf("%3d ", j);
//     }
//     printf("\n");

//     for (int i = 0; i < strlen(S1); i++) {
//         printf("%2d: ", i);
//         for (int j = 0; j < strlen(S2); j++) {
//             printf("%3d ", LCP[i][j]);
//         }
//         printf("\n");
//     }

//     // Clean up
//     for (int i = 0; i < strlen(S1); i++) {
//         free(LCP[i]);
//     }
//     free(LCP);

//     return 0;
// }