#!/usr/bin/env python3
"""
Script to download FASTQ files from md5sum.txt
Usage: python download_fastq.py
"""

import os
import sys
import hashlib
import urllib.request
from pathlib import Path


def calculate_md5(filepath):
    """Calculate MD5 checksum of a file."""
    md5_hash = hashlib.md5()
    with open(filepath, "rb") as f:
        # Read file in chunks to handle large files
        for chunk in iter(lambda: f.read(4096 * 1024), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def download_with_progress(url, output_path):
    """Download file with progress indicator."""
    def report_progress(block_num, block_size, total_size):
        downloaded = block_num * block_size
        percent = min(downloaded * 100 / total_size, 100) if total_size > 0 else 0
        mb_downloaded = downloaded / (1024 * 1024)
        mb_total = total_size / (1024 * 1024)
        print(f"\rProgress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", end="")
    
    print(f"Downloading to: {output_path}")
    urllib.request.urlretrieve(url, output_path, reporthook=report_progress)
    print()  # New line after progress


def main():
    # Get script directory and data directory
    script_dir = Path(__file__).parent
    data_dir = script_dir / "data" / "raw"
    md5sum_file = data_dir / "md5sum.txt"
    
    # Create data/raw directory if it doesn't exist
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if md5sum.txt exists
    if not md5sum_file.exists():
        print(f"Error: {md5sum_file} not found!")
        sys.exit(1)
    
    print("=" * 60)
    print("FASTQ File Download Script")
    print("=" * 60)
    print(f"Download directory: {data_dir}")
    print()
    
    # Read md5sum.txt and download files
    with open(md5sum_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header line
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        
        # Parse TSV line
        parts = line.split('\t')
        if len(parts) < 4:
            print(f"Warning: Invalid line format: {line}")
            continue
        
        filename, filesize, expected_md5, download_link = parts
        output_file = data_dir / filename
        
        # Check if file already exists
        if output_file.exists():
            print(f"File already exists: {filename}")
            print("Verifying MD5 checksum...")
            
            calculated_md5 = calculate_md5(output_file)
            
            if calculated_md5 == expected_md5:
                print("✓ MD5 checksum verified. Skipping download.")
            else:
                print("✗ MD5 checksum mismatch! Re-downloading...")
                print(f"Expected: {expected_md5}")
                print(f"Got: {calculated_md5}")
                output_file.unlink()
                download_with_progress(download_link, output_file)
                
                # Verify after re-download
                calculated_md5 = calculate_md5(output_file)
                if calculated_md5 == expected_md5:
                    print("✓ MD5 checksum verified after re-download.")
                else:
                    print("✗ MD5 checksum still mismatched!")
                    print(f"Expected: {expected_md5}")
                    print(f"Got: {calculated_md5}")
        else:
            print(f"Downloading: {filename}")
            print(f"Size: {filesize} bytes")
            
            try:
                download_with_progress(download_link, output_file)
                
                # Verify MD5 checksum after download
                print("Verifying MD5 checksum...")
                calculated_md5 = calculate_md5(output_file)
                
                if calculated_md5 == expected_md5:
                    print("✓ MD5 checksum verified.")
                else:
                    print("✗ MD5 checksum mismatch!")
                    print(f"Expected: {expected_md5}")
                    print(f"Got: {calculated_md5}")
                    print("Warning: File may be corrupted!")
            except Exception as e:
                print(f"Error downloading {filename}: {e}")
                if output_file.exists():
                    output_file.unlink()
        
        print()
    
    print("=" * 60)
    print("Download completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
