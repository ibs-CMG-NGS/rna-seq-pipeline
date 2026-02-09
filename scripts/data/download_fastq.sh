#!/bin/bash

# Script to download FASTQ files from md5sum.txt
# Usage: ./download_fastq.sh

set -e  # Exit on error

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATA_DIR="$SCRIPT_DIR/data/raw"
MD5SUM_FILE="$DATA_DIR/md5sum.txt"

# Create data/raw directory if it doesn't exist
mkdir -p "$DATA_DIR"

# Check if md5sum.txt exists
if [ ! -f "$MD5SUM_FILE" ]; then
    echo "Error: $MD5SUM_FILE not found!"
    exit 1
fi

echo "========================================"
echo "FASTQ File Download Script"
echo "========================================"
echo "Download directory: $DATA_DIR"
echo ""

# Read md5sum.txt and download files (skip header line)
tail -n +2 "$MD5SUM_FILE" | while IFS=$'\t' read -r filename filesize md5sum download_link; do
    # Skip empty lines
    if [ -z "$filename" ]; then
        continue
    fi
    
    output_file="$DATA_DIR/$filename"
    
    # Check if file already exists
    if [ -f "$output_file" ]; then
        echo "File already exists: $filename"
        echo "Verifying MD5 checksum..."
        
        # Calculate MD5 of existing file
        if command -v md5sum &> /dev/null; then
            calculated_md5=$(md5sum "$output_file" | awk '{print $1}')
        elif command -v md5 &> /dev/null; then
            calculated_md5=$(md5 -q "$output_file")
        else
            echo "Warning: md5sum/md5 command not found. Skipping verification."
            continue
        fi
        
        if [ "$calculated_md5" = "$md5sum" ]; then
            echo "✓ MD5 checksum verified. Skipping download."
        else
            echo "✗ MD5 checksum mismatch! Re-downloading..."
            rm -f "$output_file"
            wget -O "$output_file" "$download_link"
        fi
    else
        echo "Downloading: $filename"
        echo "Size: $filesize bytes"
        echo "URL: $download_link"
        
        # Download with wget
        wget -O "$output_file" "$download_link"
        
        # Verify MD5 checksum after download
        echo "Verifying MD5 checksum..."
        if command -v md5sum &> /dev/null; then
            calculated_md5=$(md5sum "$output_file" | awk '{print $1}')
        elif command -v md5 &> /dev/null; then
            calculated_md5=$(md5 -q "$output_file")
        else
            echo "Warning: md5sum/md5 command not found. Skipping verification."
            continue
        fi
        
        if [ "$calculated_md5" = "$md5sum" ]; then
            echo "✓ MD5 checksum verified."
        else
            echo "✗ MD5 checksum mismatch!"
            echo "Expected: $md5sum"
            echo "Got: $calculated_md5"
            echo "Warning: File may be corrupted!"
        fi
    fi
    echo ""
done

echo "========================================"
echo "Download completed!"
echo "========================================"
