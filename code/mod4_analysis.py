"""
Mod 4 Independence Analysis for Gilbreath's Conjecture
=======================================================
Tobias Research Project

KEY INSIGHT (Odlyzko 1993):
- For k≥1, n≥2: d_k(n) is always even
- The recurrence: d_{k+1}(n) ≡ d_k(n) + d_k(n+1) (mod 4)
- This is Pascal's triangle mod 4!
- If d_1(n)/2 values are asymptotically independent mod 2,
  then Gilbreath's conjecture follows.

We investigate this independence computationally.
"""

import numpy as np
from sympy import nextprime
from collections import Counter
import time

# ============================================================
# PART 1: Generate d_1(n) = |p_{n+1} - p_n| (prime gaps)
# and study d_1(n)/2 mod 2
# ============================================================


def get_prime_gaps(n_primes):
    """Get first n prime gaps: d_1(n) = p_{n+1} - p_n"""
    primes = []
    p = 2
    for _ in range(n_primes + 1):
        primes.append(p)
        p = nextprime(p)
    gaps = [primes[i + 1] - primes[i] for i in range(n_primes)]
    return primes[:-1], gaps


def analyze_gap_mod2(gaps):
    """
    d_1(n) = prime gap. For n≥2, d_1(n) is even.
    Let b(n) = d_1(n)/2 mod 2.
    b(n) = 0 means gap ≡ 0 (mod 4)
    b(n) = 1 means gap ≡ 2 (mod 4)

    If these are "random", each should be ~50/50.
    """
    # Skip first gap (which is 1, from 2→3)
    even_gaps = gaps[1:]  # All even gaps (from p≥3)

    b_values = [(g // 2) % 2 for g in even_gaps]

    counter = Counter(b_values)
    total = len(b_values)

    print("=" * 60)
    print("d₁(n)/2 mod 2 DISTRIBUTION")
    print("=" * 60)
    print(f"Total even gaps analyzed: {total}")
    print(f"b(n)=0 (gap ≡ 0 mod 4): {counter[0]} ({counter[0] / total * 100:.2f}%)")
    print(f"b(n)=1 (gap ≡ 2 mod 4): {counter[1]} ({counter[1] / total * 100:.2f}%)")
    print(f"Deviation from 50%: {abs(counter[0] / total - 0.5) * 100:.4f}%")

    return b_values


# ============================================================
# PART 2: Pairwise Correlation Analysis
# ============================================================


def measure_correlations(b_values, max_lag=50):
    """
    KEY TEST: Are b(n) and b(n+lag) correlated?
    If truly independent, Corr(b(n), b(n+lag)) → 0.

    We measure:
    1. Pr(b(n)=b(n+lag)) - should be 0.5 if independent
    2. Autocorrelation function
    """
    n = len(b_values)
    b = np.array(b_values, dtype=float)

    # Center the values
    mean_b = np.mean(b)
    var_b = np.var(b)
    b_centered = b - mean_b

    print("\n" + "=" * 60)
    print("PAIRWISE CORRELATION ANALYSIS")
    print("=" * 60)
    print(f"Mean of b(n): {mean_b:.6f} (should be ~0.5)")
    print(f"Variance: {var_b:.6f} (should be ~0.25)")

    print(f"\n{'Lag':>5} {'Pr(same)':>10} {'AutoCorr':>12} {'|Deviation|':>12}")
    print("-" * 45)

    correlations = []
    for lag in range(1, max_lag + 1):
        # Probability they're the same
        same = sum(1 for i in range(n - lag) if b_values[i] == b_values[i + lag])
        pr_same = same / (n - lag)

        # Autocorrelation
        if var_b > 0:
            autocorr = np.mean(b_centered[: n - lag] * b_centered[lag:]) / var_b
        else:
            autocorr = 0.0

        correlations.append((lag, pr_same, autocorr))

        if lag <= 20 or lag % 10 == 0:
            print(
                f"{lag:5d} {pr_same:10.6f} {autocorr:12.6f} {abs(pr_same - 0.5):12.6f}"
            )

    return correlations


# ============================================================
# PART 3: Higher-order independence (k-tuples)
# ============================================================


def test_tuple_independence(b_values, tuple_size=3):
    """
    Test if consecutive k-tuples of b(n) values are uniformly distributed.
    For k bits, there are 2^k possible patterns.
    If independent, each should appear with probability 1/2^k.
    """
    n = len(b_values)

    print(f"\n{'=' * 60}")
    print(f"{tuple_size}-TUPLE INDEPENDENCE TEST")
    print(f"{'=' * 60}")

    n_patterns = 2**tuple_size
    expected = (n - tuple_size + 1) / n_patterns

    # Count all k-tuples
    pattern_counts = Counter()
    for i in range(n - tuple_size + 1):
        pattern = tuple(b_values[i : i + tuple_size])
        pattern_counts[pattern] += 1

    print(f"Expected count per pattern: {expected:.1f}")
    print(f"Number of patterns: {n_patterns}")
    print()

    # Chi-squared test
    chi_sq = 0
    for pattern in sorted(pattern_counts.keys()):
        count = pattern_counts[pattern]
        deviation = (count - expected) ** 2 / expected
        chi_sq += deviation
        pct = count / (n - tuple_size + 1) * 100
        expected_pct = 100 / n_patterns
        print(f"  {pattern}: {count:7d} ({pct:.3f}%, expected {expected_pct:.3f}%)")

    # Check for missing patterns
    for i in range(n_patterns):
        pattern = tuple(int(x) for x in format(i, f"0{tuple_size}b"))
        if pattern not in pattern_counts:
            print(f"  {pattern}: MISSING!")

    print(f"\nChi-squared statistic: {chi_sq:.4f}")
    print(f"Degrees of freedom: {n_patterns - 1}")
    # Rough critical value at 5% significance
    critical = n_patterns - 1 + 2 * np.sqrt(2 * (n_patterns - 1))
    print(f"Rough critical value (5%): {critical:.1f}")
    print(
        f"Independent? {'YES (likely)' if chi_sq < critical else 'NO (correlation detected)'}"
    )

    return chi_sq, pattern_counts


# ============================================================
# PART 4: The Pascal's Triangle Mod 4 Connection
# ============================================================


def pascal_mod4_analysis(gaps, max_depth=200):
    """
    Odlyzko's key equation: d_{k+1}(n) ≡ d_k(n) + d_k(n+1) (mod 4)

    Since d_k(n) is always 0 or 2 (mod 4) for k≥1, n≥2,
    we can track the evolution as binary: d_k(n)/2 mod 2.

    The rule becomes: b_{k+1}(n) = b_k(n) XOR b_k(n+1)
    (since (0+0)/2=0, (0+2)/2=1, (2+0)/2=1, (2+2)/2=2≡0 mod 2)

    Wait - but the actual values can be > 2! The mod 4 reduction
    tells us the PARITY structure. Let's verify this computationally.
    """
    print("\n" + "=" * 60)
    print("PASCAL'S TRIANGLE MOD 4 VERIFICATION")
    print("=" * 60)

    # Start with actual d_1 values (prime gaps, all even for n≥2)
    d1 = [g for g in gaps[1:]]  # skip gap from 2→3

    # Track actual differences and their mod 4 values
    current_actual = d1[:500]  # use first 500 even gaps
    current_mod4 = [x % 4 for x in current_actual]

    print(f"\nRow 1 (first 20 d₁(n)): {current_actual[:20]}")
    print(f"Row 1 mod 4:            {current_mod4[:20]}")

    # Now iterate: d_{k+1}(n) = |d_k(n) - d_k(n+1)|
    # Check: does d_{k+1}(n) mod 4 ≡ (d_k(n) + d_k(n+1)) mod 4?

    violations = 0
    total_checks = 0

    for depth in range(min(max_depth, len(current_actual) - 1)):
        next_actual = [
            abs(current_actual[i + 1] - current_actual[i])
            for i in range(len(current_actual) - 1)
        ]
        next_mod4 = [x % 4 for x in next_actual]

        # Verify Odlyzko's congruence
        for i in range(len(next_actual)):
            predicted_mod4 = (current_mod4[i] + current_mod4[i + 1]) % 4
            actual_mod4 = next_mod4[i]
            total_checks += 1
            if predicted_mod4 != actual_mod4:
                violations += 1

        if depth < 10:
            print(f"Row {depth + 2} mod 4 (first 15): {next_mod4[:15]}")

        current_actual = next_actual
        current_mod4 = next_mod4

        if len(current_actual) < 2:
            break

    print(f"\nVerification: {total_checks} checks, {violations} violations")
    print(f"Odlyzko's congruence holds: {'YES' if violations == 0 else 'NO'}")

    return violations == 0


# ============================================================
# PART 5: XOR Structure and First Column
# ============================================================


def xor_first_column_analysis(b_values, max_depth=300):
    """
    Since d_{k+1}(n) mod 4 ≡ d_k(n) + d_k(n+1) mod 4,
    and d_k(n)/2 mod 2 follows XOR rule:

    b_{k+1}(n) = b_k(n) XOR b_k(n+1)

    But this is ONLY true mod 4. The actual values could differ.

    Let's track the XOR triangle and compare to the actual triangle.
    The first column of the XOR triangle determines d_k(1) mod 4.
    If d_k(1) ≡ 2 (mod 4) for all k, then d_k(1) is never 0 or ≥4,
    which would mean d_k(1) = 2, so d_{k-1}(1) = |d_{k-1}(1) - 2| or similar...

    Actually, what we need: d_k(1) = 1 for all k.
    Since d_k(1) is odd (always), d_k(1) mod 4 = 1 or 3.
    If mod 4 structure shows d_k(1) mod 4 is "random" (1 or 3),
    combined with small absolute values, d_k(1) = 1 follows.
    """
    print("\n" + "=" * 60)
    print("XOR TRIANGLE vs ACTUAL TRIANGLE (First Column)")
    print("=" * 60)

    # Build XOR triangle from b(n) values
    current_b = list(b_values[:max_depth])
    xor_first_col = [current_b[0]]

    for depth in range(len(current_b) - 1):
        next_b = [current_b[i] ^ current_b[i + 1] for i in range(len(current_b) - 1)]
        if len(next_b) == 0:
            break
        xor_first_col.append(next_b[0])
        current_b = next_b

    # Count values in first column
    counter = Counter(xor_first_col)
    print(f"XOR first column (first 50): {xor_first_col[:50]}")
    print(f"\nDistribution: 0s={counter[0]}, 1s={counter[1]}")
    print(f"Fraction of 1s: {counter[1] / len(xor_first_col):.4f}")

    # The XOR first column tells us d_k(1)/2 mod 2 (for the EVEN part)
    # But d_k(1) is always odd! So we need to think differently...

    # Actually: d_0(1)=2 (first prime), d_1(1)=1 (|3-2|=1), and for k≥1:
    # d_k(1) is always ODD.
    # d_k(n) for n≥2 is always EVEN (for k≥1).
    # So d_k(1) = |d_{k-1}(1) - d_{k-1}(2)| = |odd - even| = odd. ✓

    # The mod 4 structure of d_k(2), d_k(3), ... determines how d_k(1) evolves.

    return xor_first_col


# ============================================================
# PART 6: Correlation with distance — the decorrelation rate
# ============================================================


def decorrelation_study(b_values):
    """
    KEY QUESTION: How fast does correlation decay?

    If Corr(b(n), b(n+lag)) ~ C/lag^alpha for some alpha,
    that would be a new empirical result.

    Faster decay → stronger independence → closer to proving Gilbreath.
    """
    print("\n" + "=" * 60)
    print("DECORRELATION RATE STUDY")
    print("=" * 60)

    n = len(b_values)
    b = np.array(b_values, dtype=float)
    mean_b = np.mean(b)
    var_b = np.var(b)
    b_centered = b - mean_b

    lags = list(range(1, 201))
    autocorrs = []

    for lag in lags:
        autocorr = np.mean(b_centered[: n - lag] * b_centered[lag:]) / var_b
        autocorrs.append(abs(autocorr))

    # Fit power law: |autocorr(lag)| ~ C * lag^(-alpha)
    # Use log-log regression
    log_lags = np.log(np.array(lags[5:], dtype=float))  # skip first few
    log_acorr = np.log(np.array(autocorrs[5:]) + 1e-10)

    # Simple linear regression
    valid = np.isfinite(log_acorr)
    if np.sum(valid) > 10:
        coeffs = np.polyfit(log_lags[valid], log_acorr[valid], 1)
        alpha = -coeffs[0]
        C_fit = np.exp(coeffs[1])

        print(f"\nPower law fit: |Corr(lag)| ~ {C_fit:.6f} * lag^(-{alpha:.4f})")
        print(f"Decay exponent α = {alpha:.4f}")
        print()

        if alpha > 0.5:
            print("→ Moderate/fast decorrelation — encouraging for independence!")
        elif alpha > 0:
            print("→ Slow decorrelation — some residual dependence persists")
        else:
            print("→ No clear decorrelation — independence not supported")

    # Print some key correlations
    print(f"\n{'Lag':>5} {'|AutoCorr|':>12}")
    print("-" * 20)
    for i, lag in enumerate(lags):
        if lag in [1, 2, 3, 5, 10, 20, 50, 100, 150, 200]:
            print(f"{lag:5d} {autocorrs[i]:12.6f}")

    return lags, autocorrs


# ============================================================
# PART 7: Small primes dependency (Odlyzko's observation)
# ============================================================


def small_prime_dependency(primes, gaps):
    """
    Odlyzko noted: among d_1(n)/2 for 2 ≤ n ≤ 5001,
    2956 of 5000 are odd → 59.1% (excess of odds for small n).

    This was noted by Roesler as dependency among small primes.
    Let's track how this excess evolves with n.
    """
    print("\n" + "=" * 60)
    print("SMALL PRIME DEPENDENCY (Odlyzko's observation)")
    print("=" * 60)
    print("Tracking excess of odd d₁(n)/2 values as n grows")
    print()

    even_gaps = gaps[1:]  # skip first gap
    b_values = [(g // 2) % 2 for g in even_gaps]

    # Track running fraction of 1s
    checkpoints = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000]

    running_sum = 0
    print(f"{'n':>8} {'%odd (b=1)':>12} {'excess':>10} {'1/√n bound':>12}")
    print("-" * 45)

    for i, b in enumerate(b_values):
        running_sum += b
        n = i + 1

        if n in checkpoints or (n > 50000 and n % 50000 == 0):
            pct_odd = running_sum / n
            excess = pct_odd - 0.5
            bound = 1 / np.sqrt(n)
            print(f"{n:8d} {pct_odd:12.6f} {excess:+10.6f} {bound:12.6f}")

    # Final
    n = len(b_values)
    pct_odd = running_sum / n
    print(f"\nFinal ({n} values): {pct_odd:.6f} ({pct_odd * 100:.3f}% odd)")
    print(f"Expected if random: 0.500000")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("🔬 MOD 4 INDEPENDENCE ANALYSIS FOR GILBREATH'S CONJECTURE")
    print("=" * 60)

    N_PRIMES = 100000
    print(f"\nGenerating {N_PRIMES} primes...")
    t0 = time.time()
    primes, gaps = get_prime_gaps(N_PRIMES)
    print(f"Done in {time.time() - t0:.1f}s. Largest prime: {primes[-1]}")

    # Part 1: Basic distribution
    b_values = analyze_gap_mod2(gaps)

    # Part 2: Pairwise correlations
    correlations = measure_correlations(b_values, max_lag=50)

    # Part 3: k-tuple independence
    for k in [2, 3, 4, 5]:
        test_tuple_independence(b_values, tuple_size=k)

    # Part 4: Pascal's triangle mod 4
    pascal_mod4_analysis(gaps)

    # Part 5: XOR first column
    xor_first_col = xor_first_column_analysis(b_values)

    # Part 6: Decorrelation rate
    decorrelation_study(b_values)

    # Part 7: Small prime dependency
    small_prime_dependency(primes, gaps)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("=" * 60)
