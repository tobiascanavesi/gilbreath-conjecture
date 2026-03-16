#!/usr/bin/env python3
"""Quick verification that all data files are valid and complete."""

import csv
from pathlib import Path

basedir = Path(__file__).parent.parent / 'data'

files_to_check = [
    ('bias_decay_data.csv', 8, ['scale', 'N', 'epsilon', 'log_x', 'epsilon_times_logx']),
    ('autocorrelation_data.csv', 80, ['scale', 'lag', 'correlation']),
    ('xor_convergence_data.csv', 50, ['m', 'theoretical_bias_0.14', 'theoretical_bias_0.07', 'theoretical_bias_0.04']),
    ('jk_size_data.csv', 1000, ['k', 'jk_size']),
]

print("="*80)
print("DATA FILE VERIFICATION")
print("="*80)

all_good = True
for filename, expected_rows, expected_cols in files_to_check:
    filepath = basedir / filename
    print(f"\n{filename}:")
    
    if not filepath.exists():
        print(f"  ✗ FILE MISSING")
        all_good = False
        continue
    
    # Read and verify
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    # Check column headers
    if reader.fieldnames != expected_cols:
        print(f"  ✗ Column mismatch")
        print(f"    Expected: {expected_cols}")
        print(f"    Got: {reader.fieldnames}")
        all_good = False
        continue
    
    # Check row count
    if len(rows) != expected_rows:
        print(f"  ✗ Row count mismatch: expected {expected_rows}, got {len(rows)}")
        all_good = False
        continue
    
    # Check for numeric values
    sample_row = rows[0]
    all_numeric = True
    for col in expected_cols:
        try:
            if col in ['scale', 'N', 'k', 'm', 'lag']:
                int(sample_row[col])
            else:
                float(sample_row[col])
        except ValueError:
            all_numeric = False
            break
    
    if not all_numeric:
        print(f"  ✗ Non-numeric values found")
        all_good = False
        continue
    
    print(f"  ✓ {len(rows)} rows")
    print(f"  ✓ Columns: {', '.join(expected_cols)}")
    print(f"  ✓ All numeric values valid")

print("\n" + "="*80)
if all_good:
    print("✓ ALL FILES VERIFIED SUCCESSFULLY")
else:
    print("✗ SOME FILES HAVE ISSUES")
print("="*80)
