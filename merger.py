from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import os

def merge_sequences(input_file, batch_size=100):
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        return

    # Generate output filename based on input filename
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_merged.fasta"

    records = list(SeqIO.parse(input_file, "fasta"))
    merged_records = []

    for i in range(0, len(records), batch_size):
        batch = records[i:i + batch_size]
        
        merged_seq = "".join(str(record.seq) for record in batch)
        
        merged_record = SeqRecord(
            Seq(merged_seq),
            id=f"merged_{i//batch_size + 1}",
            description=f"Merged sequences {i+1} to {i+len(batch)}"
        )
        
        merged_records.append(merged_record)

    # Write all merged records to one output file
    SeqIO.write(merged_records, output_file, "fasta")
    print(f"Saved {len(merged_records)} merged sequences to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Merge FASTA sequences in batches')
    parser.add_argument('input_file', help='Path to the input FASTA file')
    parser.add_argument('--batch-size', type=int, default=100, 
                      help='Number of sequences to merge in each batch (default: 100)')
    
    args = parser.parse_args()
    merge_sequences(args.input_file, args.batch_size)

if __name__ == "__main__":
    main()