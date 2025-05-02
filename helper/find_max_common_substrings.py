from Bio import SeqIO
from itertools import combinations
import sys

def hamming_distance(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))

def find_maximal_common_substrings(s1, s2, k, tau):
    max_subs = []
    len1, len2 = len(s1), len(s2)

    for i in range(len1):
        for j in range(len2):
            mismatches = 0
            l = 0
            max_len = min(len1 - i, len2 - j)

            while l < max_len:
                if s1[i + l] != s2[j + l]:
                    mismatches += 1
                    if mismatches > k:
                        break
                l += 1

            if l < tau:
                continue

            # Maximality: check if 1-char extension stays within k mismatches
            left_mismatch = (
                i > 0 and j > 0 and
                hamming_distance(s1[i - 1:i + l], s2[j - 1:j + l]) <= k
            )
            right_mismatch = (
                i + l < len1 and j + l < len2 and
                hamming_distance(s1[i:i + l + 1], s2[j:j + l + 1]) <= k
            )

            if not left_mismatch and not right_mismatch:
                substr1 = s1[i:i + l]
                substr2 = s2[j:j + l]
                max_subs.append((i, j, l, substr1, substr2))

    return max_subs

def process_fasta(filename, k, tau):
    records = list(SeqIO.parse(filename, "fasta"))
    total_count = 0

    for idx1, idx2 in combinations(range(len(records)), 2):
        s1 = str(records[idx1].seq)
        s2 = str(records[idx2].seq)
        common_subs = find_maximal_common_substrings(s1, s2, k, tau)
        for pos1, pos2, length, substr1, substr2 in common_subs:
            print(f"Pair: {idx1} {idx2} Position1: {pos1} Position2: {pos2} Length: {length} Common Substring: {substr1}/{substr2}")
        total_count += len(common_subs)

    print(f"\nnumber of all maximal common substrings up to {k} mismatches (length ≥ {tau}): {total_count}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python find_max_common_substrings.py <input.fasta> <k> <tau>")
    else:
        fasta_file = sys.argv[1]
        k = int(sys.argv[2])
        tau = int(sys.argv[3])
        process_fasta(fasta_file, k, tau)