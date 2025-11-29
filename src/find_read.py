#!/usr/bin/env python3
import gzip

def find_problem_read(filename, target_id):
    """Find specific read in FASTQ file"""
    print(f"Searching for {target_id} in {filename}...")
    with gzip.open(filename, 'rt') as f:
        read_num = 0
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            
            read_num += 1
            
            if target_id in header:
                print(f"\nFound at read number: {read_num}")
                print(f"Header: {header.strip()}")
                print(f"Seq length: {len(seq.strip())}")
                print(f"Seq: {seq.strip()[:100]}")
                print(f"Plus: {plus.strip()}")
                print(f"Qual length: {len(qual.strip())}")
                print(f"Qual: {qual.strip()[:100]}")
                
                if len(seq.strip()) != len(qual.strip()):
                    print(f"\n⚠️  LENGTH MISMATCH!")
                else:
                    print(f"\n✓ Lengths match")
                return
            
            if read_num % 1000000 == 0:
                print(f"Scanned {read_num:,} reads...")
    
    print(f"Read {target_id} not found after scanning {read_num:,} reads")

# Search for the problematic read
find_problem_read('results/trimmed/Ctrl_3_1.fastq.gz', '1:1101:18385:1098')
