#!/usr/bin/env python3
"""
Restricted k-t Longest Common Substring (Rkt-LCS) with match reporting + matched substring.
"""

import sys

def read_fasta(path):
    """Read a FASTA file; return list of (seq_id, sequence) tuples."""
    seqs = []
    with open(path) as f:
        header = None
        buf = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    seqs.append((header, ''.join(buf)))
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if header is not None:
            seqs.append((header, ''.join(buf)))
    return seqs

def first_approx_match_info(sub, seq, k):
    """
    Return (position, matched_substring) where `matched_substring` from `seq`
    has Hamming distance <= k to `sub`. If none, return (None, None).
    """
    L = len(sub)
    for i in range(len(seq) - L + 1):
        mismatches = 0
        for a, b in zip(sub, seq[i:i+L]):
            if a != b:
                mismatches += 1
                if mismatches > k:
                    break
        if mismatches <= k:
            return i, seq[i:i+L]
    return None, None

def restricted_k_t_lcs_report(fasta_path, k, t, tau):
    sequences = read_fasta(fasta_path)

    for seq_id, seq in sequences:
        best = None  
        # best: (start, length, substring, [ (other_id, pos, matched_sub), ... ] )
        n = len(seq)

        # try lengths from n down to tau
        for L in range(n, tau-1, -1):
            found = False
            for start in range(0, n - L + 1):
                u = seq[start:start+L]
                matches = []
                # look for up to t other-sequence matches
                for other_id, other_seq in sequences:
                    if other_id == seq_id:
                        continue
                    pos, matched = first_approx_match_info(u, other_seq, k)
                    if pos is not None:
                        matches.append((other_id, pos, matched))
                        if len(matches) >= t:
                            break
                if len(matches) >= t:
                    best = (start, L, u, matches[:t])
                    found = True
                    break
            if found:
                break

        # print results
        print(f"Sequence {seq_id} ({seq_id}):")
        if best:
            start, L, u, matches = best
            print(f"  Position: {start}")
            print(f"  Length: {L}")
            print(f"  Substring: {u}")
            print("  Matches in other sequences:")
            for other_id, pos, matched in matches:
                print(f"    - Sequence {other_id}: position {pos} ({matched})")
        else:
            print("  No valid substring found")
        print("-" * 40)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python rkt_lcs_report.py <FASTA> <k> <t> <tau>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    k = int(sys.argv[2])
    t = int(sys.argv[3])
    tau = int(sys.argv[4])

    restricted_k_t_lcs_report(fasta_file, k, t, tau)