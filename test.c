#include "dependencies/include/mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "flouri_lcp_table.h" // Include table.h


void test_LCP(const char* s1, const char* s2, int k) {
    int l1 = strlen(s1), l2 = strlen(s2);
    printf("Comparing LCP s1 = \"%s\" with s2 = \"%s\" k = \"%d\"\n", s1, s2, k);
    int** lcpTable = compute_k_LCP(s1, s2, k);
    print2DArray(lcpTable, l1, l2);
    free2DArray(lcpTable, l1);
}

void test_LengthStat(char* S[], int m, int k, int i, int p) {
    int** LCP_i = (int**) malloc((m - 1) * sizeof(int*));
    for (int j = i; j < m; j++) {
        LCP_i[j - i] = compute_k_LCP_max(S[i], S[j], k);
        printf("\n=== LCP_i[%d] ===\n", j);
        printArray(LCP_i[j - i], strlen(S[j]));
    }
    
    // Call compute_LengthStats.
    int** LengthStat = compute_LengthStat(LCP_i, S, m, p, i);
    int L = strlen(S[i]) - 1 - p; 
    printf("\n=== LengthStat Table ===\n");
    print2DArray(LengthStat, L+1, m+1);
}

void check_result(const char *testName, char *result, const char *expected) {
    if (result == NULL && expected == NULL) {
        printf("%s PASSED: result = NULL, expected = NULL\n", testName);
    } else if (result != NULL && strcmp(result, expected) == 0) {
        printf("%s PASSED: result = \"%s\", expected = \"%s\"\n", testName, result, expected);
    } else {
        printf("%s FAILED: result = \"%s\", expected = \"%s\"\n", testName, result, expected);
    }
    free(result);  /* Assume Rkt_LCS allocated memory with malloc/calloc */
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

    char* S1[] = {"abc", "abd", "ccc"};
    int m = 3;          // m+1 = 3 strings in S
    int i_fixed = 0;    // Fixed index i (we use S[0] as the reference)
    int p = 1;          // Choose p (must be less than strlen(S[i_fixed]))
    int k = 2; 
    // test_LengthStat(S1, m, k, i_fixed, p);
    // test_LengthStat(S2, 3, 1, 0, 2);

    // const char* S2[] = {"TTGAC", "CGAAAT", "TGGTA"};
    // Rkt_LCS(S2, 3, 1, 3);

    /* --------------------------- */
    /* Test 2: Using S2            */
    /* --------------------------- */
    {
        char* S2[] = {"TTGAC", "CGAAAT", "TGGTA"};
        int m = 3;    /* Number of strings */
        int k = 1;    /* Parameter k (interpretation depends on your implementation) */
        int t = 3;    /* Parameter t */
        
        char *result = Rkt_LCS(S2, m, k, t);
        char *expected = "TGA";  
        
        check_result("Test 2", result, expected);
    }

    /* --------------------------- */
    /* Test 3: A different example */
    /* --------------------------- */
    {
        char* S3[] = {"ABCDEF", "ZBCXYZ", "WBCPQR", "BCSTUV"};
        int m = 4;
        int k = 1;
        int t = 4;
        
        char *result = Rkt_LCS(S3, m, k, t);
        char *expected = "BCD";  

        // test_LengthStat(S3, m, k, 0, 1);
        
        check_result("Test 3", result, expected);
    }
    
    /* --------------------------- */
    /* Test 4: All strings identical */
    /* --------------------------- */
    {
        char* S4[] = {"HELLO", "HELLO", "HELLO"};
        int m = 3;
        int k = 1;
        int t = 3;
        
        /* If all strings are identical, the LCS is expected to be the full string */
        char *result = Rkt_LCS(S4, m, k, t);
        char *expected = "HELLO";
        
        check_result("Test 4", result, expected);
    }

    /* --------------------------- */
    /* Test 5: Can't find common string */
    /* --------------------------- */
    {
        char* S5[] = {"HELLO", "HELLO", "AA"};
        int m = 3;
        int k = 0;
        int t = 3;
        
        char *result = Rkt_LCS(S5, m, k, t);
        char *expected = NULL;
        
        check_result("Test 5", result, expected);
    }

    /* --------------------------- */
    /* Test 6*/
    /* --------------------------- */
    {
        char* S6[] = {"aaaaa", "bbbbb", "ccccc", "ddcdd"};
        int m = 4;
        int k = 2;
        int t = 2;
        
        char *result = Rkt_LCS(S6, m, k, t);
        char *expected = "ccc";
        check_result("Test 6", result, expected);
    }

    return 0;
}