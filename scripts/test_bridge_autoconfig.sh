#!/bin/bash
# Test bridge script with auto-config

echo "============================================================"
echo "Test 1: Bridge with auto-config (no --config flag)"
echo "============================================================"

cd /data_3tb/shared/rna-seq-pipeline

# Remove existing config to test auto-generation
rm -f config/projects/paths_mouse_chd8.yaml

# Run bridge with just project-id
conda run -n rna-seq-pipeline python scripts/bridge_to_de_pipeline.py \
    --project-id mouse-chd8 \
    --skip-de \
    --yes

echo ""
echo "============================================================"
echo "Test 2: Bridge with existing config (should reuse)"
echo "============================================================"

# Run again - should use existing config
conda run -n rna-seq-pipeline python scripts/bridge_to_de_pipeline.py \
    --project-id mouse-chd8 \
    --skip-de \
    --yes

echo ""
echo "============================================================"
echo "Test 3: Bridge with explicit config"
echo "============================================================"

# Run with explicit config
conda run -n rna-seq-pipeline python scripts/bridge_to_de_pipeline.py \
    --config config/projects/paths_mouse_chd8.yaml \
    --project-id mouse-chd8 \
    --skip-de \
    --yes
