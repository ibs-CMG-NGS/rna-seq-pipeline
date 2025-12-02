#!/usr/bin/env python3
import gzip
import sys

def check_fastq(filename, max_reads=1000):
    """Check FASTQ file format"""
    print(f"Checking {filename}...")
    with gzip.open(filename, 'rt') as f:
        for i in range(max_reads):
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            
            if not header.startswith('@'):
                print(f"ERROR at read {i+1}: Header doesn't start with @")
                print(f"Header: {header[:80]}")
                return False
            
            if not plus.startswith('+'):
                print(f"ERROR at read {i+1}: Plus line doesn't start with +")
                print(f"Plus line: {plus[:80]}")
                print(f"Previous header: {header[:80]}")
                return False
            
            seq_len = len(seq.strip())
            qual_len = len(qual.strip())
            if seq_len != qual_len:
                print(f"ERROR at read {i+1}: Length mismatch")
                print(f"Header: {header.strip()}")
                print(f"Seq length: {seq_len}, Qual length: {qual_len}")
                print(f"Seq: {seq[:80]}")
                print(f"Qual: {qual[:80]}")
                return False
    
    print(f"✓ Checked {i+1} reads - all OK!")
    return True

# Check both R1 and R2
check_fastq('results/trimmed/Ctrl_3_1.fastq.gz', max_reads=50000)
print()
check_fastq('results/trimmed/Ctrl_3_2.fastq.gz', max_reads=50000)
