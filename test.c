#include "dependencies/include/mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "flouri_lcp_table.h" // Include table.h


//---------------------------------------------------------
// Helper function: print2DArray
//
// Prints a 2D integer array with the specified number of rows and columns.
// Uses C's printf function.
//
void print2DArray(int** arr, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        printf("Row %d: ", i);
        for (int j = 0; j < cols; j++) {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }
}



//---------------------------------------------------------
// Helper function: free2DArray
//
// Frees a 2D array that was allocated using new[]. 
//
void free2DArray(int** arr, int rows) {
    for (int i = 0; i < rows; i++) {
        free(arr[i]);
    }
    free(arr);
}


void test_LCP(const char* s1, const char* s2, int k) {
    int l1 = strlen(s1), l2 = strlen(s2);
    printf("Comparing LCP s1 = \"%s\" with s2 = \"%s\" k = \"%d\"\n", s1, s2, k);
    int** lcpTable = compute_k_LCP(s1, s2, k);
    print2DArray(lcpTable, l1, l2);
    free2DArray(lcpTable, l1);
}

void test_LengthStat(char* S[], int m, int k, int i, int p) {
    int*** LCP_i = (int***) malloc((m - 1) * sizeof(int**));
    for (int j = 0; j < m; j++) {
        LCP_i[j] = compute_k_LCP(S[i], S[j], k);
    }
    
    // Call compute_LengthStats.
    int** LengthStat = compute_LengthStats(LCP_i, S, m, p, i);
    int L = strlen(S[i]) - 1 - p; 
    printf("\n=== LengthStat Table ===\n");
    print2DArray(LengthStat, L+1, m+1);
}


//---------------------------------------------------------
// Main test routine
//
int main() {
    // --- Test Case 1: Test and print the result of compute_k_LCP ---
    char* s1 = "abc";
    char* s2 = "abd";
    char* s3 = "GTACAAT";
    char* s4 = "CTTGTA";
    // test_LCP(s1, s2, 1);
    // test_LCP(s3, s4, 2);

    const char* S1[] = {"abc", "abd", "ccc"};
    int m = 3;          // m+1 = 3 strings in S
    int i_fixed = 0;    // Fixed index i (we use S[0] as the reference)
    int p = 1;          // Choose p (must be less than strlen(S[i_fixed]))
    int k = 2; // For example, k = 3 if S[0] is "abc"
    // test_LengthStat(S1, m, k, i_fixed, p);
    const char* S2[] = {"TTGAC", "CGAAAT", "TGGTA"};
    // test_LengthStat(S2, 3, 1, 0, 2);
    Rkt_LCS(S2, 3, 1, 3);

    return 0;
}