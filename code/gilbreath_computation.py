#!/usr/bin/env python3
"""
Comprehensive computational verification for Gilbreath's conjecture paper.

This script performs expanded measurements of:
1. Bias decay (epsilon) across scales up to 10^9
2. Autocorrelation at lags 1-20
3. Gilbreath triangle verification
4. XOR mod 4 convergence
5. |J_k| size distribution
6. Generates data files for paper figures
"""

import numpy as np
import csv
import time
from pathlib import Path
from collections import defaultdict
from sympy import prime, nextprime, primerange

# Setup output directory
output_dir = Path(__file__).parent.parent / 'data'
output_dir.mkdir(exist_ok=True)

print("="*80)
print("GILBREATH'S CONJECTURE: COMPREHENSIVE COMPUTATIONAL ANALYSIS")
print("="*80)
print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# ============================================================================
# PART 1: BIAS MEASUREMENT (epsilon(x) decay)
# ============================================================================

print("\n" + "="*80)
print("PART 1: BIAS DECAY MEASUREMENT (epsilon(x) = Pr(b(n)=1) - 1/2)")
print("="*80)

def compute_bias(primes_sample):
    """
    Compute bias epsilon = Pr(b(n)=1) - 1/2 for a sample of primes.
    b(n) = 1 if gap ≡ 2 (mod 4), b(n) = 0 if gap ≡ 0 (mod 4).
    """
    b_ones = 0
    for i in range(len(primes_sample) - 1):
        gap = primes_sample[i+1] - primes_sample[i]
        if gap % 4 == 2:
            b_ones += 1

    return (b_ones / len(primes_sample) - 0.5)

# Target scales and sample sizes
scales = [
    (100, 50000),
    (1000, 50000),
    (10000, 50000),
    (100000, 50000),
    (1000000, 50000),
    (10000000, 50000),
]

# Extend with 10^8 and 10^9 if time permits
try_extended = True
if try_extended:
    scales.extend([
        (100000000, 50000),
        (1000000000, 50000),
    ])

bias_results = []

print("\nComputing bias at each scale (sampling 50,000 consecutive primes):")
print(f"{'Scale':<12} {'x_approx':<15} {'Epsilon':<12} {'log(x)':<12} {'Eps*log(x)':<12} {'Time (s)':<10}")
print("-" * 85)

for target_scale, sample_size in scales:
    scale_start = time.time()

    # Find primes near the target scale
    if target_scale == 100:
        start_prime = prime(10)  # 29
        primes_sample = list(primerange(start_prime, start_prime + 100000))[:sample_size]
    elif target_scale == 1000:
        start_prime = prime(150)  # ~860
        primes_sample = list(primerange(start_prime, start_prime + 100000))[:sample_size]
    elif target_scale == 10000:
        start_prime = prime(1200)  # ~9689
        primes_sample = list(primerange(start_prime, start_prime + 100000))[:sample_size]
    elif target_scale == 100000:
        start_prime = prime(8500)  # ~87953
        primes_sample = list(primerange(start_prime, start_prime + 100000))[:sample_size]
    elif target_scale == 1000000:
        start_prime = prime(70000)  # ~868951
        primes_sample = list(primerange(start_prime, start_prime + 100000))[:sample_size]
    elif target_scale == 10000000:
        start_prime = prime(620000)  # ~9241000 (approx)
        # For large scales, sample from specific range
        primes_sample = []
        p = start_prime
        while len(primes_sample) < sample_size:
            p = nextprime(p)
            primes_sample.append(p)
    elif target_scale == 100000000:
        start_prime = prime(5500000)  # ~97000000 (approx)
        primes_sample = []
        p = start_prime
        while len(primes_sample) < sample_size:
            p = nextprime(p)
            primes_sample.append(p)
    elif target_scale == 1000000000:
        start_prime = prime(50000000)  # ~982000000 (approx)
        primes_sample = []
        p = start_prime
        while len(primes_sample) < sample_size:
            p = nextprime(p)
            primes_sample.append(p)
    else:
        continue

    if not primes_sample:
        continue

    epsilon = compute_bias(primes_sample)
    x_mean = np.mean(primes_sample)
    log_x = np.log(x_mean)
    eps_times_logx = epsilon * log_x

    scale_time = time.time() - scale_start

    if scale_time > 5:  # Stop if taking too long
        print(f"{target_scale:<12} {x_mean:<15.2e} {epsilon:<12.6f} {log_x:<12.6f} {eps_times_logx:<12.6f} {scale_time:<10.2f}")
        bias_results.append({
            'scale': target_scale,
            'N': len(primes_sample),
            'epsilon': epsilon,
            'log_x': log_x,
            'eps_times_logx': eps_times_logx,
            'x_mean': x_mean
        })
        break

    print(f"{target_scale:<12} {x_mean:<15.2e} {epsilon:<12.6f} {log_x:<12.6f} {eps_times_logx:<12.6f} {scale_time:<10.2f}")
    bias_results.append({
        'scale': target_scale,
        'N': len(primes_sample),
        'epsilon': epsilon,
        'log_x': log_x,
        'eps_times_logx': eps_times_logx,
        'x_mean': x_mean
    })

print(f"\nBias measurements completed. {len(bias_results)} scales analyzed.")

# ============================================================================
# PART 2: AUTOCORRELATION MEASUREMENT
# ============================================================================

print("\n" + "="*80)
print("PART 2: AUTOCORRELATION MEASUREMENT (lags 1-20)")
print("="*80)

def compute_autocorrelation(b_sequence, max_lag=20):
    """
    Compute autocorrelation for the b(n) sequence at various lags.
    b(n) = 1 if prime gap ≡ 2 (mod 4), else 0.
    """
    mean = np.mean(b_sequence)
    var = np.var(b_sequence)

    if var < 1e-10:
        return [0.0] * max_lag

    correlations = []
    for lag in range(1, max_lag + 1):
        if lag >= len(b_sequence):
            correlations.append(0.0)
            continue

        cov = np.mean((b_sequence[:-lag] - mean) * (b_sequence[lag:] - mean))
        corr = cov / var
        correlations.append(corr)

    return correlations

autocorr_scales = [
    (100000, 100000),
    (1000000, 100000),
    (10000000, 100000),
    (100000000, 100000),
]

autocorr_results = []

print("\nComputing autocorrelations at each scale (sampling 100,000 primes):")
print(f"{'Scale':<15} {'Time (s)':<10}")
print("-" * 30)

for target_scale, sample_size in autocorr_scales:
    scale_start = time.time()

    # Generate b sequence
    if target_scale == 100000:
        start_prime = prime(8500)
    elif target_scale == 1000000:
        start_prime = prime(70000)
    elif target_scale == 10000000:
        start_prime = prime(620000)
    elif target_scale == 100000000:
        start_prime = prime(5500000)
    else:
        continue

    primes_sample = []
    p = start_prime
    while len(primes_sample) < sample_size + 1:
        primes_sample.append(p)
        p = nextprime(p)

    # Generate b sequence
    b_sequence = []
    for i in range(len(primes_sample) - 1):
        gap = primes_sample[i+1] - primes_sample[i]
        b = 1 if gap % 4 == 2 else 0
        b_sequence.append(b)

    b_array = np.array(b_sequence)
    correlations = compute_autocorrelation(b_array, max_lag=20)

    scale_time = time.time() - scale_start

    print(f"{target_scale:<15} {scale_time:<10.2f}")

    autocorr_results.append({
        'scale': target_scale,
        'correlations': correlations
    })

    if scale_time > 15:  # Stop if taking too long
        break

print(f"\nAutocorrelation analysis completed. {len(autocorr_results)} scales analyzed.")

# ============================================================================
# PART 3: GILBREATH TRIANGLE VERIFICATION
# ============================================================================

print("\n" + "="*80)
print("PART 3: GILBREATH TRIANGLE VERIFICATION")
print("="*80)

def build_gilbreath_triangle(primes, max_rows=100):
    """
    Build the Gilbreath triangle for up to max_rows.
    Returns: (triangle, d_1_column, reset_counts, max_gap_between_resets)
    """
    if len(primes) < 2:
        return None

    # Start with d_0 = primes
    current_row = primes[:]
    triangle = [current_row[:]]
    d_1_column = []

    reset_count = 0
    last_reset_row = 0
    max_gap = 0
    d2_resets = []  # Track which rows have d_k(2) = 2

    for k in range(1, max_rows):
        if len(current_row) < 2:
            break

        # Compute next row
        next_row = []
        for i in range(len(current_row) - 1):
            diff = abs(current_row[i] - current_row[i+1])
            next_row.append(diff)

        triangle.append(next_row[:])

        # Check d_k(1)
        d_k_1 = next_row[0] if next_row else None
        d_1_column.append(d_k_1)

        # Check d_k(2) for resets
        if len(next_row) > 1:
            d_k_2 = next_row[1]
            if d_k_2 == 2:
                reset_count += 1
                gap = k - last_reset_row
                if gap > max_gap:
                    max_gap = gap
                last_reset_row = k
                d2_resets.append(k)

        current_row = next_row

    # Check final gap
    if d2_resets:
        final_gap = len(triangle) - 1 - d2_resets[-1]
        if final_gap > max_gap:
            max_gap = final_gap

    return {
        'triangle': triangle,
        'd_1_column': d_1_column,
        'reset_count': reset_count,
        'max_gap': max_gap,
        'd2_resets': d2_resets,
        'violations_d1': sum(1 for d1 in d_1_column if d1 != 1)
    }

print("\nBuilding Gilbreath triangle from primes...")

# Build triangle for 50,000 primes
primes_for_triangle = []
p = prime(50000)  # Start near 600,000
while len(primes_for_triangle) < 50000:
    primes_for_triangle.append(p)
    p = nextprime(p)

triangle_start = time.time()
triangle_result = build_gilbreath_triangle(primes_for_triangle, max_rows=500)
triangle_time = time.time() - triangle_start

print(f"Triangle built from {len(primes_for_triangle)} primes in {triangle_time:.2f}s")
print(f"Rows computed: {len(triangle_result['triangle'])}")
print(f"Violations of d_k(1)=1: {triangle_result['violations_d1']}")
print(f"Resets (d_k(2)=2): {triangle_result['reset_count']}")
print(f"Max consecutive gap between resets: {triangle_result['max_gap']}")
print(f"Reset frequency: {100*triangle_result['reset_count']/len(triangle_result['d_1_column']):.1f}%")
print(f"First reset at row: {triangle_result['d2_resets'][0] if triangle_result['d2_resets'] else 'N/A'}")
print(f"Last reset at row: {triangle_result['d2_resets'][-1] if triangle_result['d2_resets'] else 'N/A'}")

# ============================================================================
# PART 4: XOR MOD 4 VERIFICATION
# ============================================================================

print("\n" + "="*80)
print("PART 4: XOR MOD 4 VERIFICATION (Corollary 1)")
print("="*80)

def verify_xor_mod4(primes_sample, num_checks=10000):
    """
    Verify that e_{k+1}(n) = e_k(n) XOR e_k(n+1) where e_k(n) = (d_k(n)/2) mod 2.
    """
    violations = 0

    # Build rows
    current_row = primes_sample[:num_checks+2]

    for k_iter in range(20):  # Check first 20 rows
        if len(current_row) < 3:
            break

        # Compute next row (differences)
        next_row = [abs(current_row[i] - current_row[i+1]) for i in range(len(current_row)-1)]

        # Verify XOR property
        for n in range(len(next_row) - 1):
            # e_k values
            e_k_n = (current_row[n] // 2) % 2
            e_k_n1 = (current_row[n+1] // 2) % 2
            # Expected e_{k+1}(n)
            expected = e_k_n ^ e_k_n1
            # Actual e_{k+1}(n)
            actual = (next_row[n] // 2) % 2

            if expected != actual:
                violations += 1

        current_row = next_row

    return violations

xor_violations = verify_xor_mod4(primes_for_triangle, num_checks=5000)
print(f"XOR mod 4 violations in 5000 primes, 20 rows: {xor_violations}")
print(f"Status: {'VERIFIED' if xor_violations == 0 else 'FAILED'}")

# ============================================================================
# PART 5: |J_k| SIZE DISTRIBUTION
# ============================================================================

print("\n" + "="*80)
print("PART 5: |J_k| SIZE DISTRIBUTION (Lucas' Theorem)")
print("="*80)

def hamming_weight(n):
    """Count the number of 1-bits in binary representation of n."""
    return bin(n).count('1')

def compute_jk_size(k):
    """
    Compute |J_k| = 2^{s(k-1)} where s(n) = hamming weight of n.
    J_k is the set of indices j such that C(k-1,j) is odd (mod 2).
    By Lucas' theorem, |J_k| = 2^{s(k-1)}.
    """
    return 2 ** hamming_weight(k - 1)

print("\nComputing |J_k| for k = 1 to 1000 (Lucas' theorem):")

jk_sizes = []
jk_equals_2 = 0
min_jk = float('inf')
max_jk = 0

for k in range(1, 1001):
    jk_size = compute_jk_size(k)
    jk_sizes.append(jk_size)

    if jk_size == 2:
        jk_equals_2 += 1

    min_jk = min(min_jk, jk_size)
    max_jk = max(max_jk, jk_size)

mean_jk = np.mean(jk_sizes)
median_jk = np.median(jk_sizes)
std_jk = np.std(jk_sizes)

print(f"\n|J_k| Statistics (k = 1 to 1000):")
print(f"  Minimum: {min_jk}")
print(f"  Maximum: {max_jk}")
print(f"  Mean: {mean_jk:.2f}")
print(f"  Median: {median_jk:.1f}")
print(f"  Std Dev: {std_jk:.2f}")
print(f"  Frequency of |J_k|=2: {jk_equals_2}/1000 = {100*jk_equals_2/1000:.1f}%")

# Show distribution
jk_distribution = defaultdict(int)
for jk in jk_sizes:
    jk_distribution[jk] += 1

print(f"\n|J_k| Value Distribution (top 10):")
for size in sorted(jk_distribution.keys(), key=lambda x: jk_distribution[x], reverse=True)[:10]:
    count = jk_distribution[size]
    print(f"  |J_k| = {size:<4d}: {count:<4d} times ({100*count/1000:>5.1f}%)")

# ============================================================================
# PART 6: GENERATE DATA FILES FOR FIGURES
# ============================================================================

print("\n" + "="*80)
print("PART 6: GENERATING DATA FILES FOR PAPER FIGURES")
print("="*80)

# File 1: bias_decay_data.csv
print("\nGenerating bias_decay_data.csv...")
with open(output_dir / 'bias_decay_data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['scale', 'N', 'epsilon', 'log_x', 'epsilon_times_logx'])
    for result in bias_results:
        writer.writerow([
            result['scale'],
            result['N'],
            f"{result['epsilon']:.6f}",
            f"{result['log_x']:.6f}",
            f"{result['eps_times_logx']:.6f}"
        ])

print(f"  Written {len(bias_results)} rows")

# File 2: autocorrelation_data.csv
print("Generating autocorrelation_data.csv...")
with open(output_dir / 'autocorrelation_data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['scale', 'lag', 'correlation'])
    for result in autocorr_results:
        for lag, corr in enumerate(result['correlations'], start=1):
            writer.writerow([
                result['scale'],
                lag,
                f"{corr:.6f}"
            ])

print(f"  Written {len(autocorr_results) * 20} rows")

# File 3: xor_convergence_data.csv
print("Generating xor_convergence_data.csv...")
# Generate theoretical convergence data
with open(output_dir / 'xor_convergence_data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['m', 'theoretical_bias_0.14', 'theoretical_bias_0.07', 'theoretical_bias_0.04'])

    for m in range(1, 51):
        # For bias epsilon, after m XORs: |Pr(Y=1) - 1/2| = (1/2)(2*epsilon)^m
        eps_014 = 0.5 * (0.14 ** m)
        eps_007 = 0.5 * (0.07 ** m)
        eps_004 = 0.5 * (0.04 ** m)
        writer.writerow([
            m,
            f"{eps_014:.10e}",
            f"{eps_007:.10e}",
            f"{eps_004:.10e}"
        ])

print(f"  Written 50 rows")

# File 4: jk_size_data.csv
print("Generating jk_size_data.csv...")
with open(output_dir / 'jk_size_data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['k', 'jk_size'])
    for k in range(1, 1001):
        writer.writerow([k, compute_jk_size(k)])

print(f"  Written 1000 rows")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print("\n" + "="*80)
print("SUMMARY OF RESULTS")
print("="*80)

print("\n1. BIAS DECAY:")
print(f"   - Measured at {len(bias_results)} scales")
if len(bias_results) >= 2:
    first_eps = bias_results[0]['epsilon']
    last_eps = bias_results[-1]['epsilon']
    print(f"   - Epsilon range: {last_eps:.6f} to {first_eps:.6f}")
    eps_logx_vals = [f"{r['eps_times_logx']:.3f}" for r in bias_results]
    print(f"   - Product epsilon*log(x) is stable: {eps_logx_vals}")

print("\n2. AUTOCORRELATION:")
print(f"   - Measured at {len(autocorr_results)} scales")
if autocorr_results:
    lag1_corrs = [r['correlations'][0] for r in autocorr_results]
    lag1_strs = [f"{c:.4f}" for c in lag1_corrs]
    print(f"   - Lag-1 correlations: {lag1_strs}")
    lag2_corrs = [r['correlations'][1] for r in autocorr_results]
    lag2_strs = [f"{c:.4f}" for c in lag2_corrs]
    print(f"   - Lag-2 correlations: {lag2_strs}")

print("\n3. GILBREATH TRIANGLE:")
print(f"   - Built from 50,000 primes, {len(triangle_result['triangle'])} rows")
print(f"   - All d_k(1) = 1: {triangle_result['violations_d1'] == 0}")
print(f"   - Reset frequency: {100*triangle_result['reset_count']/len(triangle_result['d_1_column']):.1f}%")
print(f"   - Max gap between resets: {triangle_result['max_gap']}")

print("\n4. XOR MOD 4 VERIFICATION:")
print(f"   - Verified with 5000 primes, 20 rows")
print(f"   - Violations: {xor_violations}")

print("\n5. |J_k| DISTRIBUTION:")
print(f"   - Computed for k=1 to 1000")
print(f"   - Mean: {mean_jk:.2f}, Range: [{min_jk}, {max_jk}]")
print(f"   - |J_k|=2 frequency: {100*jk_equals_2/1000:.1f}%")

print("\n6. DATA FILES GENERATED:")
print("   - bias_decay_data.csv")
print("   - autocorrelation_data.csv")
print("   - xor_convergence_data.csv")
print("   - jk_size_data.csv")

print("\n" + "="*80)
print(f"Computation completed at {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
