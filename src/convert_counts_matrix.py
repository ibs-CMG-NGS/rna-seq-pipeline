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
import pandas as pd
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
    
    # Read the featureCounts output, skipping comment lines
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep='\t', comment='#')
    
    # Display original columns
    print(f"\nOriginal columns: {list(df.columns)}")
    print(f"Original shape: {df.shape}")
    
    # Keep only Geneid and sample count columns
    # Remove: Chr, Start, End, Strand, Length
    metadata_cols = ['Chr', 'Start', 'End', 'Strand', 'Length']
    cols_to_keep = ['Geneid'] + [col for col in df.columns if col not in metadata_cols + ['Geneid']]
    
    clean_df = df[cols_to_keep].copy()
    
    # Clean up sample names (remove path and .bam extension)
    # e.g., "results/aligned/GABA_8/Aligned.sortedByCoord.out.bam" -> "GABA_8"
    new_columns = ['Geneid']
    for col in clean_df.columns[1:]:
        if '/' in col:
            # Extract sample name from path
            sample_name = col.split('/')[-2] if col.endswith('.bam') else col.split('/')[-1]
            sample_name = sample_name.replace('Aligned.sortedByCoord.out.bam', '').strip()
        else:
            sample_name = col.replace('.bam', '')
        new_columns.append(sample_name)
    
    clean_df.columns = new_columns
    
    # Set Geneid as index
    clean_df.set_index('Geneid', inplace=True)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save to CSV
    print(f"\nSaving cleaned matrix to: {output_file}")
    clean_df.to_csv(output_file)
    
    print(f"\nCleaned count matrix shape: {clean_df.shape}")
    print(f"Sample names: {list(clean_df.columns)}")
    print("\nFirst few rows:")
    print(clean_df.head())
    
    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Total genes: {len(clean_df)}")
    print(f"Total samples: {len(clean_df.columns)}")
    print(f"\nGenes with zero counts across all samples: {(clean_df.sum(axis=1) == 0).sum()}")
    print(f"Average counts per gene: {clean_df.sum(axis=1).mean():.2f}")
    
    return clean_df


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
