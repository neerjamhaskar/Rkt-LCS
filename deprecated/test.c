#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define MAX_SEQ 1000000   // max number of sequences
#define STR_LEN 51        // fixed length for all sequences
#define MAX_LINE 1024

// Compute Hamming distance
int hamming_distance(const char* a, const char* b, int len) {
    if (!a || !b) {
        fprintf(stderr, "Error: NULL pointer in hamming_distance\n");
        return -1;
    }
    int dist = 0;
    for (int i = 0; i < len; i++) {
        if (a[i] != b[i]) dist++;
    }
    return dist;
}

// Print substring to stream
void print_substring(FILE* out, const char* s, int start, int len) {
    if (!out || !s) {
        fprintf(stderr, "Error: NULL pointer in print_substring\n");
        return;
    }
    for (int i = 0; i < len && s[start + i]; i++) {
        fputc(s[start + i], out);
    }
}

int read_fasta(const char* filename, char sequences[][STR_LEN + 1]) {
    if (!filename || !sequences) {
        fprintf(stderr, "Error: NULL pointer in read_fasta\n");
        return -1;
    }

    printf("=== Starting read_fasta ===\n");
    printf("Attempting to open file: %s\n", filename);
    
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file %s: %s\n", filename, strerror(errno));
        return -1;
    }

    printf("Successfully opened file %s\n", filename);
    char line[MAX_LINE];
    int count = 0;
    int line_num = 0;

    // Initialize all sequences to empty strings
    for (int i = 0; i < MAX_SEQ; i++) {
        sequences[i][0] = '\0';
    }

    while (fgets(line, sizeof(line), file)) {
        line_num++;
        
        // Skip empty lines
        if (line[0] == '\n' || line[0] == '\r') {
            continue;
        }
        
        // If this is a header line, skip it
        if (line[0] == '>' || line[0] == '1' || line[0] == '2' || line[0] == '3') {
            printf("Found header at line %d: %s", line_num, line);
            continue;
        }

        // This must be a sequence line
        size_t len = strlen(line);
        if (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[len - 1] = '\0';
            len--;
        }
        // Handle Windows-style line endings
        if (len > 0 && line[len - 1] == '\r') {
            line[len - 1] = '\0';
            len--;
        }

        printf("Reading sequence at line %d (length=%zu): %s\n", line_num, len, line);

        if (len != STR_LEN) {
            fprintf(stderr, "Invalid sequence length at line %d (%zu, expected %d): \"%s\"\n", 
                    line_num, len, STR_LEN, line);
            fclose(file);
            return -1;
        }

        if (count >= MAX_SEQ) {
            fprintf(stderr, "Too many sequences (>%d) at line %d\n", MAX_SEQ, line_num);
            fclose(file);
            return -1;
        }

        strncpy(sequences[count], line, STR_LEN);
        sequences[count][STR_LEN] = '\0';  // Ensure null termination
        printf("Stored sequence %d: %s\n", count + 1, sequences[count]);
        count++;
    }

    printf("Finished reading file. Found %d sequences.\n", count);
    fclose(file);
    return count;
}

// Check if a substring match is maximal (can't be extended while maintaining k mismatches)
int is_maximal_match(const char* s1, const char* s2, int start1, int start2, int len, int k) {
    // Check if we can extend to the left
    if (start1 > 0 && start2 > 0) {
        int left_dist = 0;
        if (s1[start1 - 1] != s2[start2 - 1]) left_dist = 1;
        
        // Count mismatches in current window
        int curr_dist = 0;
        for (int i = 0; i < len && s1[start1 + i] && s2[start2 + i]; i++) {
            if (s1[start1 + i] != s2[start2 + i]) curr_dist++;
        }
        
        // If extending left maintains exactly k mismatches, it's not maximal
        if (curr_dist + left_dist == k) return 0;
    }
    
    // Check if we can extend to the right
    if (s1[start1 + len] && s2[start2 + len]) {
        int right_dist = 0;
        if (s1[start1 + len] != s2[start2 + len]) right_dist = 1;
        
        // Count mismatches in current window
        int curr_dist = 0;
        for (int i = 0; i < len && s1[start1 + i] && s2[start2 + i]; i++) {
            if (s1[start1 + i] != s2[start2 + i]) curr_dist++;
        }
        
        // If extending right maintains exactly k mismatches, it's not maximal
        if (curr_dist + right_dist == k) return 0;
    }
    
    return 1;  // Cannot extend in either direction while maintaining k mismatches
}

void find_all_k_mismatches(char sequences[][STR_LEN + 1], int m, int k, int min_len, const char* outfile) {
    if (!sequences || !outfile || m <= 0 || k < 0 || min_len <= 0) {
        fprintf(stderr, "Error: Invalid arguments to find_all_k_mismatches\n");
        return;
    }

    printf("Entering find_all_k_mismatches with m=%d, k=%d, min_len=%d\n", m, k, min_len);
    printf("Finding substrings with exactly %d mismatches\n", k);
    
    FILE* fout = fopen(outfile, "w");
    if (!fout) {
        fprintf(stderr, "Error opening output file %s: %s\n", outfile, strerror(errno));
        return;
    }
    printf("Successfully opened output file: %s\n", outfile);

    // First pass: find the maximum length among all matches with exactly k mismatches
    int max_match_length = 0;
    for (int i = 0; i < m; i++) {
        if (!sequences[i][0]) continue;
        for (int j = i + 1; j < m; j++) {
            if (!sequences[j][0]) continue;
            for (int len = min_len; len <= STR_LEN; len++) {
                for (int pi = 0; pi <= STR_LEN - len; pi++) {
                    for (int pj = 0; pj <= STR_LEN - len; pj++) {
                        int dist = hamming_distance(sequences[i] + pi, sequences[j] + pj, len);
                        if (dist == k && len > max_match_length) {
                            max_match_length = len;
                        }
                    }
                }
            }
        }
    }

    // Reset file position to start writing matches
    rewind(fout);
    
    int total_matches = 0;
    int maximal_matches = 0;
    int longest_maximal_matches = 0;
    
    fprintf(fout, "Substring pairs with exactly %d mismatches:\n", k);
    fprintf(fout, "(M = maximal, cannot be extended while maintaining %d mismatches)\n", k);
    fprintf(fout, "(L = longest maximal, maximal and of maximum length %d)\n\n", max_match_length);
    
    for (int i = 0; i < m; i++) {
        if (!sequences[i][0]) {
            fprintf(stderr, "Error: Empty sequence at index %d\n", i);
            continue;
        }
        
        for (int j = i + 1; j < m; j++) {
            if (!sequences[j][0]) {
                fprintf(stderr, "Error: Empty sequence at index %d\n", j);
                continue;
            }
            
            for (int len = min_len; len <= STR_LEN; len++) {
                for (int pi = 0; pi <= STR_LEN - len; pi++) {
                    for (int pj = 0; pj <= STR_LEN - len; pj++) {
                        int dist = hamming_distance(sequences[i] + pi, sequences[j] + pj, len);
                        if (dist == -1) {
                            fprintf(stderr, "Error computing Hamming distance\n");
                            continue;
                        }
                        if (dist == k) {
                            int is_maximal = is_maximal_match(sequences[i], sequences[j], pi, pj, len, k);
                            int is_longest = (len == max_match_length);
                            
                            fprintf(fout, "s%d s%d %d %d %d ", i + 1, j + 1, pi, pj, len);
                            print_substring(fout, sequences[i], pi, len);
                            fprintf(fout, " ");
                            print_substring(fout, sequences[j], pj, len);
                            
                            if (is_maximal) {
                                if (is_longest) {
                                    fprintf(fout, " (L)");
                                    longest_maximal_matches++;
                                } else {
                                    fprintf(fout, " (M)");
                                }
                                maximal_matches++;
                            }
                            fprintf(fout, "\n");
                            total_matches++;
                        }
                    }
                }
            }
        }
    }

    fprintf(fout, "\nSummary:\n");
    fprintf(fout, "Total matching pairs with exactly %d mismatches: %d\n", k, total_matches);
    fprintf(fout, "Number of maximal matching pairs (cannot be extended): %d\n", maximal_matches);
    fprintf(fout, "Number of longest maximal matching pairs (length %d): %d\n", max_match_length, longest_maximal_matches);
    
    printf("\nSummary:\n");
    printf("Total matching pairs with exactly %d mismatches: %d\n", k, total_matches);
    printf("Number of maximal matching pairs (cannot be extended): %d\n", maximal_matches);
    printf("Number of longest maximal matching pairs (length %d): %d\n", max_match_length, longest_maximal_matches);
    
    fclose(fout);
}

int main(int argc, char *argv[]) {
    printf("=== Program Start ===\n");
    
    if (!argv) {
        fprintf(stderr, "Error: NULL argv pointer\n");
        return 1;
    }
    
    printf("Command line arguments received: %d\n", argc - 1);
    for (int i = 0; i < argc; i++) {
        if (!argv[i]) {
            fprintf(stderr, "Error: NULL argument at index %d\n", i);
            return 1;
        }
        printf("argv[%d] = %s\n", i, argv[i]);
    }
    
    if (argc != 4) {
        fprintf(stderr, "\nError: Wrong number of arguments\n");
        fprintf(stderr, "Usage: %s <fasta_file> <min_len> <k>\n", argv[0]);
        fprintf(stderr, "  <fasta_file>: Path to the FASTA file containing sequences\n");
        fprintf(stderr, "  <min_len>: Minimum length of matching substrings (1-%d)\n", STR_LEN);
        fprintf(stderr, "  <k>: Maximum number of allowed mismatches\n");
        fprintf(stderr, "\nExample: %s test.txt 15 1\n", argv[0]);
        return 1;
    }

    char (*sequences)[STR_LEN + 1] = calloc(MAX_SEQ, sizeof(char[STR_LEN + 1]));
    if (!sequences) {
        fprintf(stderr, "Error: Failed to allocate memory for sequences: %s\n", strerror(errno));
        return 1;
    }

    const char* fasta_file = argv[1];
    int min_len = atoi(argv[2]);
    int k = atoi(argv[3]);
    const char* output_file = "output.txt";

    printf("\nInput parameters:\n");
    printf("  fasta_file: %s\n", fasta_file);
    printf("  min_len: %d\n", min_len);
    printf("  k: %d\n", k);
    printf("  output_file: %s\n\n", output_file);

    // Validate parameters
    if (min_len <= 0 || min_len > STR_LEN) {
        fprintf(stderr, "Error: min_len must be between 1 and %d\n", STR_LEN);
        free(sequences);
        return 1;
    }
    if (k < 0) {
        fprintf(stderr, "Error: k must be non-negative\n");
        free(sequences);
        return 1;
    }

    printf("Opening file: %s\n", fasta_file);
    int num_seqs = read_fasta(fasta_file, sequences);
    if (num_seqs < 0) {
        fprintf(stderr, "Failed to read sequences from %s.\n", fasta_file);
        free(sequences);
        return 1;
    }
    if (num_seqs < 2) {
        fprintf(stderr, "Need at least 2 sequences to find matches.\n");
        free(sequences);
        return 1;
    }

    printf("Read %d sequences from %s\n", num_seqs, fasta_file);
    printf("Parameters: min_len=%d, k=%d\n", min_len, k);
    
    // Validate all sequences
    for (int i = 0; i < num_seqs; i++) {
        size_t len = strlen(sequences[i]);
        if (len != STR_LEN) {
            fprintf(stderr, "Error: Sequence %d has invalid length %zu (expected %d)\n", 
                    i+1, len, STR_LEN);
            free(sequences);
            return 1;
        }
        printf("Sequence %d: %s (length: %zu)\n", i+1, sequences[i], len);
    }

    find_all_k_mismatches(sequences, num_seqs, k, min_len, output_file);

    free(sequences);
    return 0;
}