#!/usr/bin/env python3
"""
Convert featureCounts output to standard count matrix format for DE analysis.

This script removes metadata columns (Chr, Start, End, Strand, Length) 
and keeps only gene IDs and sample counts, making it compatible with 
DESeq2, edgeR, and limma-voom workflows.

Usage:
    python convert_counts_matrix.py <input_file> <output_file>
    
Example:
    python src/convert_counts_matrix.py results/counts/counts_matrix.txt results/counts/counts_matrix_clean.csv
"""

import sys
import os


def convert_featurecounts_to_clean_matrix(input_file, output_file):
    """
    Convert featureCounts output to clean count matrix.
    
    Parameters:
    -----------
    input_file : str
        Path to featureCounts output file (counts_matrix.txt)
    output_file : str
        Path to save cleaned count matrix
    """
    
    print(f"Reading input file: {input_file}")
    
    # Read the file line by line
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Skip comment lines and find header
    data_lines = []
    header_line = None
    
    for line in lines:
        if line.startswith('#'):
            continue
        if header_line is None:
            header_line = line.strip()
        else:
            data_lines.append(line.strip())
    
    # Parse header
    header = header_line.split('\t')
    print(f"\nOriginal columns: {header[:10]}...")  # Show first 10 columns
    print(f"Total columns: {len(header)}")
    
    # Find column indices to keep (Geneid and sample columns)
    # Skip: Chr, Start, End, Strand, Length (columns 1-5)
    metadata_cols = {'Chr', 'Start', 'End', 'Strand', 'Length'}
    
    keep_indices = []
    new_header = []
    
    for i, col_name in enumerate(header):
        if col_name == 'Geneid':
            keep_indices.append(i)
            new_header.append('Geneid')
        elif col_name not in metadata_cols:
            keep_indices.append(i)
            # Clean up sample name
            sample_name = col_name
            if '/' in sample_name:
                # Extract sample name from path
                # e.g., "results/aligned/GABA_8/Aligned.sortedByCoord.out.bam" -> "GABA_8"
                parts = sample_name.split('/')
                sample_name = parts[-2] if sample_name.endswith('.bam') else parts[-1]
            sample_name = sample_name.replace('Aligned.sortedByCoord.out.bam', '').replace('.bam', '').strip()
            new_header.append(sample_name)
    
    print(f"\nNew column names: {new_header}")
    print(f"Number of samples: {len(new_header) - 1}")
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Write cleaned matrix to CSV
    print(f"\nSaving cleaned matrix to: {output_file}")
    
    with open(output_file, 'w') as out:
        # Write header
        out.write(','.join(new_header) + '\n')
        
        # Write data
        for line in data_lines:
            fields = line.split('\t')
            selected_fields = [fields[i] for i in keep_indices]
            out.write(','.join(selected_fields) + '\n')
    
    # Calculate simple statistics
    total_genes = len(data_lines)
    print(f"\n=== Summary Statistics ===")
    print(f"Total genes: {total_genes}")
    print(f"Total samples: {len(new_header) - 1}")
    
    # Show first few lines
    print("\nFirst 5 genes:")
    for i, line in enumerate(data_lines[:5]):
        fields = line.split('\t')
        gene_id = fields[0]
        counts = [fields[j] for j in keep_indices[1:6]]  # Show first 5 sample counts
        print(f"  {gene_id}: {', '.join(counts)}...")
    
    print(f"\n✓ Conversion completed successfully!")
    
    return total_genes


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    try:
        convert_featurecounts_to_clean_matrix(input_file, output_file)
        print("\n✓ Conversion completed successfully!")
        print(f"\nYou can now use '{output_file}' for DE analysis with DESeq2/edgeR/limma-voom.")
        
    except Exception as e:
        print(f"\n✗ Error during conversion: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
