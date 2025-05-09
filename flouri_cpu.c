#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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



// Computes a MaxLCP_k array between S1 and S2,
// allowing up to k mismatches. An extra parameter tau is used so that
// starting positions in S1 with less than tau characters remaining are discarded.
// Moreover, the first tau characters are compared normally (mismatches counted)
// and only if they are processed (i.e. p reaches tau) do we continue.
// Thus, the returned length will include those tau characters.
int* compute_k_LCP_max(const char* S1, const char* S2, int k, int tau, int r) {
    int n = r;
    int m = strlen(S2);

    // Validate input parameters
    if (n < tau || m < tau) {
        fprintf(stderr, "Error: String lengths (S1=%d, S2=%d) are less than minimum match length (tau=%d)\n", 
                n, m, tau);
        return NULL;
    }

    // Only consider S1 positions that have at least tau characters remaining.
    int resultSize = n - tau + 1;
    if (resultSize <= 0) {
        fprintf(stderr, "Error: Invalid result size %d (n=%d, tau=%d)\n", resultSize, n, tau);
        return NULL;
    }

    int* LCP = (int*)malloc(resultSize * sizeof(int));
    if (!LCP) {
        fprintf(stderr, "Error: Memory allocation failed for LCP array (size=%d)\n", resultSize);
        return NULL;
    }

    // Create a queue to track mismatch positions (capacity = k)
    Queue* Q = createQueue(k);
    if (!Q) {
        fprintf(stderr, "Error: Failed to create queue\n");
        free(LCP);
        return NULL;
    }

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


// Structure to store position-length pair as key
typedef struct {
    int position;
    int length;
} PosLenKey;

// Structure to store position-length pair and its count
typedef struct {
    PosLenKey key;
    int count;
} PosLenCount;

// compute_k_LCP_max_multi
// Computes a MaxLCP array between S1 and multiple S2 strings,
// allowing up to k mismatches. For each position in S1, tracks how many S2 strings
// have a common substring of each length.
PosLenCount* compute_k_LCP_max_multi(const char* S1, char** S2_array, int num_S2, int k, int tau, int r, int* result_size) {
    int n = r;
    int max_possible_length = n - tau + 1;
    
    // Validate input parameters
    if (n < tau) {
        fprintf(stderr, "Error: S1 length (%d) is less than minimum match length (tau=%d)\n", n, tau);
        return NULL;
    }
    
    // Allocate array to store position-length counts
    PosLenCount* counts = NULL;
    *result_size = 0;
    int capacity = 0;
    
    // For each S2 string
    for (int s2_idx = 0; s2_idx < num_S2; s2_idx++) {
        const char* S2 = S2_array[s2_idx];
        
        // Update progress
        int progress = (s2_idx * 100) / num_S2;
        printf("Computing LCP results (%d%%)\r", progress);
        fflush(stdout);
        
        // Validate S2 length
        int m = strlen(S2);
        if (m < tau) {
            fprintf(stderr, "Error: S2[%d] length (%d) is less than minimum match length (tau=%d)\n", 
                    s2_idx, m, tau);
            continue;
        }
        
        // Compute LCP max array for this S2 string
        int* LCP_max = compute_k_LCP_max(S1, S2, k, tau, r);
        if (!LCP_max) {
            fprintf(stderr, "Error: Failed to compute LCP max array for S2[%d]\n", s2_idx);
            continue;
        }
        
        // For each starting position in S1 that has at least tau characters
        for (int start1 = 0; start1 < max_possible_length; start1++) {
            int length = LCP_max[start1];
            // For all lengths from tau up to the maximum found
            for (int l = tau; l <= length; l++) {
                // Check if we need to resize the array
                if (*result_size >= capacity) {
                    capacity = capacity == 0 ? 16 : capacity * 2;
                    counts = realloc(counts, capacity * sizeof(PosLenCount));
                    if (!counts) {
                        fprintf(stderr, "Error: Memory allocation failed for counts array\n");
                        free(LCP_max);
                        return NULL;
                    }
                }
                // Create the key for this position-length pair
                PosLenKey key = {start1, l};
                // Check if this position-length pair already exists
                int found = 0;
                for (int i = 0; i < *result_size; i++) {
                    if (counts[i].key.position == key.position && 
                        counts[i].key.length == key.length) {
                        counts[i].count++;
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
        
        // Free the MaxLCP array for this S2 string
        free(LCP_max);
    }
    
    if (*result_size == 0) {
        fprintf(stderr, "Warning: No matches found with the given parameters\n");
    }
    
    return counts;
}

// Rkt_LCS_single
//
// For a single string S1 and an array of S2 strings, this function finds
// the Rkt_LCS of S1 (candidate) that appears in at least t S2 strings
// (allowing up to k mismatches). It returns a PosLenKey containing the
// position and length of the longest such substring.
PosLenKey Rkt_LCS_single(const char* S1, char** S2_array, int num_S2, int k, int t, int tau, int r) {
    int result_size;
    PosLenCount* results = compute_k_LCP_max_multi(S1, S2_array, num_S2, k, tau, r, &result_size);
    
    if (!results) {
        fprintf(stderr, "Failed to compute multi-string LCP results\n");
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

// Structure to store a sequence from FASTA
typedef struct {
    char* name;
    char* sequence;
} FastaSequence;

// Function to read sequences from FASTA file
FastaSequence* read_fasta(const char* filename, int* num_sequences) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return NULL;
    }

    // First pass: count sequences and get lengths
    *num_sequences = 0;
    int current_length = 0;
    char line[1024];
    int in_sequence = 0;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') {
            (*num_sequences)++;
            in_sequence = 1;
        } else if (in_sequence) {
            current_length += strlen(line) - 1;  // -1 for newline
        }
    }

    // Allocate array for sequences
    FastaSequence* sequences = (FastaSequence*)malloc(*num_sequences * sizeof(FastaSequence));
    if (!sequences) {
        fprintf(stderr, "Memory allocation failed for sequences\n");
        fclose(file);
        return NULL;
    }

    // Second pass: read sequences
    rewind(file);
    int seq_idx = -1;
    int name_length = 0;
    int seq_length = 0;
    char* current_seq = NULL;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') {
            // Save previous sequence if exists
            if (seq_idx >= 0) {
                sequences[seq_idx].sequence = current_seq;
                current_seq[seq_length] = '\0';
            }

            // Start new sequence
            seq_idx++;
            name_length = strlen(line) - 1;  // -1 for newline
            sequences[seq_idx].name = (char*)malloc(name_length);
            strncpy(sequences[seq_idx].name, line + 1, name_length - 1);
            sequences[seq_idx].name[name_length - 1] = '\0';

            // Reset sequence length and allocate new sequence
            seq_length = 0;
            current_seq = (char*)malloc(1024);  // Initial size
            if (!current_seq) {
                fprintf(stderr, "Memory allocation failed for sequence\n");
                fclose(file);
                return NULL;
            }
        } else {
            // Append to current sequence
            int line_length = strlen(line) - 1;  // -1 for newline
            if (seq_length + line_length >= 1024) {
                // Reallocate if needed
                current_seq = realloc(current_seq, seq_length + line_length + 1);
                if (!current_seq) {
                    fprintf(stderr, "Memory reallocation failed for sequence\n");
                    fclose(file);
                    return NULL;
                }
            }
            strncpy(current_seq + seq_length, line, line_length);
            seq_length += line_length;
        }
    }

    // Save last sequence
    if (seq_idx >= 0) {
        sequences[seq_idx].sequence = current_seq;
        current_seq[seq_length] = '\0';
    }

    fclose(file);
    return sequences;
}

//Function to free FASTA sequences
void free_fasta(FastaSequence* sequences, int num_sequences) {
    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i].name);
        free(sequences[i].sequence);
    }
    free(sequences);
}

//Test main function
#ifndef SKIP_MAIN
int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 6) {
        printf("Usage: %s <fasta_file> <s1_seq_num> <k> <t> <tau>\n", argv[0]);
        printf("  fasta_file: Input FASTA file containing sequences\n");
        printf("  s1_seq_num: Sequence number to use as S1 (1-based index)\n");
        printf("  k: Maximum number of mismatches allowed\n");
        printf("  t: Minimum number of S2 strings that must contain the substring\n");
        printf("  tau: Minimum match length\n");
        return 1;
    }

    // Parse command line arguments
    const char* fasta_file = argv[1];
    int s1_seq_num = atoi(argv[2]);
    int k = atoi(argv[3]);
    int t = atoi(argv[4]);
    int tau = atoi(argv[5]);

    if (s1_seq_num < 1) {
        printf("Error: Sequence number must be positive\n");
        return 1;
    }
    if (k < 0 || t < 1 || tau < 1) {
        printf("Error: k must be non-negative, t and tau must be positive\n");
        return 1;
    }

    // Read sequences from FASTA file
    int num_sequences;
    FastaSequence* sequences = read_fasta(fasta_file, &num_sequences);
    if (!sequences) {
        return 1;
    }

    if (num_sequences < 2) {
        printf("Error: FASTA file must contain at least 2 sequences\n");
        free_fasta(sequences, num_sequences);
        return 1;
    }

    if (s1_seq_num > num_sequences) {
        printf("Error: Sequence number %d is larger than the number of sequences (%d)\n", 
               s1_seq_num, num_sequences);
        free_fasta(sequences, num_sequences);
        return 1;
    }

    // Use specified sequence as S1, rest as S2
    const char* S1 = sequences[s1_seq_num - 1].sequence;
    char** S2_array = (char**)malloc((num_sequences - 1) * sizeof(char*));
    if (!S2_array) {
        printf("Error: Memory allocation failed for S2 array\n");
        free_fasta(sequences, num_sequences);
        return 1;
    }

    // Fill S2_array with all sequences except the selected S1
    int s2_idx = 0;
    for (int i = 0; i < num_sequences; i++) {
        if (i != s1_seq_num - 1) {
            S2_array[s2_idx++] = sequences[i].sequence;
        }
    }
    int num_S2 = num_sequences - 1;

    int r = strlen(S1);  // Range parameter

    printf("Testing with parameters:\n");
    printf("S1: %s (sequence %d from %s)\n", sequences[s1_seq_num - 1].name, s1_seq_num, fasta_file);
    printf("k: %d (maximum mismatches)\n", k);
    printf("t: %d (minimum number of strings)\n", t);
    printf("tau: %d (minimum match length)\n", tau);
    printf("r: %d (range parameter)\n\n", r);

    // Now compute the multi-string LCP results
    printf("Computing LCP results (0%%)\r");
    fflush(stdout);
    
    int result_size;
    PosLenCount* results = compute_k_LCP_max_multi(S1, S2_array, num_S2, k, tau, r, &result_size);
    
    // Print the results
    printf("\nLCP Results:\n");
    printf("Total matches found: %d\n", result_size);
    printf("Position\tLength\tCount\n");
    printf("--------------------------------\n");
    for (int i = 0; i < result_size; i++) {
        printf("%8d\t%6d\t%5d\n", 
               results[i].key.position,
               results[i].key.length,
               results[i].count);
    }
    printf("--------------------------------\n\n");
    
    if (!results) {
        printf("Failed to compute multi-string LCP results\n");
        free(S2_array);
        free_fasta(sequences, num_sequences);
        return 1;
    }

    printf("Computing LCP results (100%%)\n");

    // Test Rkt_LCS_single
    printf("Finding longest common substring...\n");
    PosLenKey lcs_key = Rkt_LCS_single(S1, S2_array, num_S2, k, t, tau, r);
    
    if (lcs_key.position != -1) {
        printf("Longest common substring found at position %d with length %d\n", 
               lcs_key.position, lcs_key.length);
        printf("Substring: \"");
        for (int i = 0; i < lcs_key.length; i++) {
            printf("%c", S1[lcs_key.position + i]);
        }
        printf("\"\n");
    } else {
        printf("No common substring found that appears in at least %d strings\n", t);
    }

    // Clean up
    free(results);
    free(S2_array);
    free_fasta(sequences, num_sequences);

    return 0;
}
#endif