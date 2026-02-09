#!/usr/bin/env python3
import glob

samples = sorted(set([f.replace('_1.fastq.gz', '').split('/')[-1] 
                     for f in glob.glob('data/raw/*_1.fastq.gz')]))

print(f'Detected samples: {len(samples)}')
for s in samples:
    print(f'  - {s}')
