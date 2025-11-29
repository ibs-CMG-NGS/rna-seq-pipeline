import gzip
import sys
import os
import shutil

def fix_fastq(input_file, output_file):
    """Fix FASTQ files where sequences and qualities are split across 2 lines each"""
    read_count = 0
    with gzip.open(input_file, 'rt') as fin, gzip.open(output_file, 'wt') as fout:
        while True:
            # Read header
            header = fin.readline()
            if not header:
                break
            
            # Read sequence (2 lines)
            seq1 = fin.readline().strip()
            seq2 = fin.readline().strip()
            
            # Read plus line
            plus = fin.readline()
            
            # Read quality (2 lines)
            qual1 = fin.readline().strip()
            qual2 = fin.readline().strip()
            
            # Write in proper 4-line format
            fout.write(header)
            fout.write(seq1 + seq2 + '\n')
            fout.write(plus)
            fout.write(qual1 + qual2 + '\n')
            
            read_count += 1
            if read_count % 1000000 == 0:
                print(f"Processed {read_count:,} reads from {input_file}")
    
    print(f"Completed {input_file}: {read_count:,} total reads")

# Create backup directory if not exists
if not os.path.exists('backup'):
    os.makedirs('backup')
    print("Created backup directory")

# Fix Ctrl_3 files
print("\nBacking up and fixing Ctrl_3_1.fastq.gz...")
if not os.path.exists('backup/Ctrl_3_1.fastq.gz'):
    shutil.copy2('Ctrl_3_1.fastq.gz', 'backup/Ctrl_3_1.fastq.gz')
    print("Backup created: backup/Ctrl_3_1.fastq.gz")
fix_fastq('backup/Ctrl_3_1.fastq.gz', 'Ctrl_3_1.fastq.gz')

print("\nBacking up and fixing Ctrl_3_2.fastq.gz...")
if not os.path.exists('backup/Ctrl_3_2.fastq.gz'):
    shutil.copy2('Ctrl_3_2.fastq.gz', 'backup/Ctrl_3_2.fastq.gz')
    print("Backup created: backup/Ctrl_3_2.fastq.gz")
fix_fastq('backup/Ctrl_3_2.fastq.gz', 'Ctrl_3_2.fastq.gz')

print("\nAll files fixed successfully!")
