/*
 * k-t LCS under Hamming distance via Kraken CLI calls with k' mismatches
 *
 * Usage: ./program <k_prime> <t> <ref_fasta> <S_fasta> <db_name>
 *
 * 1. Builds a Kraken2 database from S_fasta
 * 2. For each candidate k (via binary search):
 *    a. Generates all k-length substrings of s_r and all variants with up to
 *       k_prime mismatches into a temporary FASTA
 *    b. Runs kraken2 classification on that FASTA against the database
 *    c. Parses kraken2 output to count how many distinct sequences in S contain
 *       at least one variant k-mer
 * 3. Reports the maximum k satisfying ≥ t matches
 *
 * Requires: kraken2 installed and in PATH
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_CMD 1024
#define MAX_LINE 4096

static int k_prime;  // maximum allowed mismatches
static int t;        // number of sequences to match

// Write all variants of length k of seq (with up to k_prime mismatches) to FASTA
static void generate_variants(const char *kmer, char *buf,
                              int pos, int mismatches,
                              int k, int id, FILE *fp) {
    if (mismatches > k_prime) return;
    if (pos == k) {
        fprintf(fp, ">kmer_%d_%s\n%s\n", id, buf, buf);
        return;
    }
    char orig = kmer[pos];
    // no mismatch at this position
    buf[pos] = orig;
    generate_variants(kmer, buf, pos+1, mismatches, k, id, fp);
    // try mismatches
    const char alphabet[] = "ACGT";
    for (int ai = 0; ai < 4; ++ai) {
        char c = alphabet[ai];
        if (c == orig) continue;
        buf[pos] = c;
        generate_variants(kmer, buf, pos+1, mismatches+1, k, id, fp);
    }
    buf[pos] = orig;
}

// Extract all k-mers of length k and their mismatched variants into FASTA
static void extract_and_variants(const char *seq, int k,
                                 const char *out_fasta) {
    FILE *fp = fopen(out_fasta, "w");
    if (!fp) { perror("fopen"); exit(EXIT_FAILURE); }
    int len = strlen(seq);
    char *buf = malloc(k+1);
    buf[k] = '\0';
    for (int i = 0; i + k <= len; ++i) {
        char *kmer = strndup(seq + i, k);
        // generate variants for this kmer
        generate_variants(kmer, buf, 0, 0, k, i, fp);
        free(kmer);
    }
    free(buf);
    fclose(fp);
}

// Build Kraken2 DB from S_fasta
static void build_db(const char *db_name, const char *S_fasta) {
    char cmd[MAX_CMD];
    snprintf(cmd, MAX_CMD,
             "kraken2-build --db %s --add-to-library %s && kraken2-build --db %s --build",
             db_name, S_fasta, db_name);
    if (system(cmd) != 0) {
        fprintf(stderr, "Error: failed to build Kraken2 database\n");
        exit(EXIT_FAILURE);
    }
}

// Run kraken2 on query_fasta, output to report_file
static void run_kraken(const char *db_name,
                       const char *query_fasta,
                       const char *out_report) {
    char cmd[MAX_CMD];
    snprintf(cmd, MAX_CMD,
             "kraken2 --db %s --report %s --output /dev/null %s",
             db_name, out_report, query_fasta);
    if (system(cmd) != 0) {
        fprintf(stderr, "Error: kraken2 classification failed\n");
        exit(EXIT_FAILURE);
    }
}

// Count number of distinct S sequences matching at least one variant k-mer
static int count_matches(const char *report_file) {
    FILE *fp = fopen(report_file, "r");
    if (!fp) { perror("fopen"); exit(EXIT_FAILURE); }
    char line[MAX_LINE];
    int count = 0;
    while (fgets(line, MAX_LINE, fp)) {
        int pct, num;
        // kraken2 report: pct\treads_covered\t...\tname
        if (sscanf(line, "%d%%\t%d", &pct, &num) == 2) {
            if (num > 0) count++;
        }
    }
    fclose(fp);
    return count;
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr,
                "Usage: %s <k_prime> <t> <ref_fasta> <S_fasta> <db_name>\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    k_prime    = atoi(argv[1]);
    t          = atoi(argv[2]);
    const char *ref_fasta = argv[3];
    const char *S_fasta   = argv[4];
    const char *db_name    = argv[5];

    // Read reference sequence
    FILE *rf = fopen(ref_fasta, "r");
    if (!rf) { perror("fopen"); return 1; }
    char *s_r = calloc(1,1);
    char *line = NULL;
    size_t cap = 0;
    while (getline(&line, &cap, rf) != -1) {
        if (line[0] == '>') continue;
        size_t len = strcspn(line, "\r\n");
        s_r = realloc(s_r, strlen(s_r) + len + 1);
        strncat(s_r, line, len);
    }
    free(line); fclose(rf);

    // Build Kraken2 database
    build_db(db_name, S_fasta);

    int n = strlen(s_r);
    int lo = 1, hi = n, best_k = 0;
    const char *kmers_file   = "/tmp/kraken_kmers.fa";
    const char *report_file  = "/tmp/kraken_report.txt";

    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        extract_and_variants(s_r, mid, kmers_file);
        run_kraken(db_name, kmers_file, report_file);
        int matches = count_matches(report_file);
        if (matches >= t) {
            best_k = mid;
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    if (best_k > 0) {
        printf("Best length = %d\n", best_k);
    } else {
        printf("No substring of s_r matches >= %d sequences within %d mismatches.\n", t, k_prime);
    }

    free(s_r);
    return EXIT_SUCCESS;
}
