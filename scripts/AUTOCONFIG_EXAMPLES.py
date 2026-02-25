#!/usr/bin/env python3
"""
Auto-config Usage Examples

Run these on the server where actual data exists.
"""

# Example 1: Auto-generate config from project ID (simplest)
from scripts.utils.auto_config import ensure_bridge_config

result = ensure_bridge_config("mouse-chd8")
print(result)

# Example 2: Generate config from project summary path
from scripts.utils.auto_config import generate_bridge_config
from pathlib import Path

config_path = generate_bridge_config(
    Path("/data_3tb/shared/output/mouse-chd8/project_summary.json")
)
print(f"Generated: {config_path}")

# Example 3: Force regenerate config
result = ensure_bridge_config("mouse-chd8", force=True)
print(result)

# Example 4: Command-line usage
"""
# Auto-discover from project ID
python scripts/utils/auto_config.py --project-id mouse-chd8

# Generate from specific summary file
python scripts/utils/auto_config.py /data_3tb/shared/output/mouse-chd8/project_summary.json

# Force overwrite
python scripts/utils/auto_config.py --project-id mouse-chd8 --force

# Custom DE pipeline location
python scripts/utils/auto_config.py --project-id mouse-chd8 --de-pipeline /path/to/de-pipeline
"""
