#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Queue structure definition
typedef struct {
    int* elements;
    int capacity;
    int size;
    int front;
    int rear;
} Queue;

// Queue operations
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

// Main K-LCF function
void K_LCF(const char* S1, const char* S2, int k) {
    int n = strlen(S1);
    int m = strlen(S2);
    int l = 0, r1 = 0, r2 = 0;

    // Create queue with capacity k
    Queue* Q = createQueue(k);

    for (int d = -m + 1; d <= n - 1; d++) {
        int i = ((-d > 0) ? -d : 0) + d;
        int j = (-d > 0) ? -d : 0;

        // Clear queue for new diagonal
        while (Q->size > 0) {
            dequeue(Q);
        }

        int s = 0, p = 0;

        while (p <= ((n - i < m - j) ? n - i : m - j) - 1) {
            if (S1[i + p] != S2[j + p]) {
                if (Q->size == k) {
                    s = minQueue(Q) + 1;
                    dequeue(Q);
                }
                enqueue(Q, p);
            }
            p++;

            if (p - s > l) {
                l = p - s;
                r1 = i + s;
                r2 = j + s;
            }
        }
    }

    // Clean up
    free(Q->elements);
    free(Q);

    // Result can be accessed through l (length), r1 (position in S1), and r2 (position in S2)
    printf("Length: %d\nPosition in S1: %d\nPosition in S2: %d\n", l, r1, r2);
}

int main() {
    // Example usage
    const char* S1 = "ABCDEF";
    const char* S2 = "ACDEFG";
    int k = 1;

    K_LCF(S1, S2, k);
    return 0;
}