"""
Gilbreath's Conjecture Explorer
================================
Tobias Research Toolkit

The conjecture: When you take the sequence of primes, compute successive
absolute differences, the first element of each row is ALWAYS 1.

We explore WHY this might be true by analyzing patterns in the
difference triangle.
"""

from sympy import prime
import time
from collections import Counter

# ============================================================
# PART 1: Build the Gilbreath Triangle
# ============================================================


def generate_gilbreath_triangle(n_primes=100):
    """Generate the Gilbreath difference triangle for the first n primes."""
    primes = [int(prime(i)) for i in range(1, n_primes + 1)]  # first n primes

    triangle = [primes]
    current = primes

    for _ in range(n_primes - 1):
        next_row = [abs(current[i + 1] - current[i]) for i in range(len(current) - 1)]
        triangle.append(next_row)
        current = next_row
        if len(current) <= 1:
            break

    return triangle


def print_triangle(triangle, max_rows=20, max_cols=20):
    """Pretty-print the Gilbreath triangle."""
    for i, row in enumerate(triangle[:max_rows]):
        display = row[:max_cols]
        label = f"Row {i:3d}: "
        first = f"[{display[0]}]" if display else "[]"
        rest = " ".join(str(x) for x in display[1:])
        print(f"{label}{first} {rest}{'  ...' if len(row) > max_cols else ''}")


def verify_conjecture(n_primes):
    """Verify that the first element of every row is 1 (after row 0)."""
    primes = [int(prime(i)) for i in range(1, n_primes + 1)]
    current = primes
    first_elements = [primes[0]]  # 2

    for i in range(n_primes - 1):
        next_row = [abs(current[j + 1] - current[j]) for j in range(len(current) - 1)]
        first_elements.append(next_row[0])
        current = next_row
        if len(current) <= 1:
            break

    # Row 0 starts with 2 (the first prime), all others should be 1
    all_ones = all(x == 1 for x in first_elements[1:])
    return all_ones, first_elements


# ============================================================
# PART 2: Statistical Analysis of Rows
# ============================================================


def analyze_row_distributions(triangle, max_rows=50):
    """Analyze the distribution of values in each row."""
    print("\n" + "=" * 60)
    print("DISTRIBUTION ANALYSIS OF GILBREATH ROWS")
    print("=" * 60)

    for i, row in enumerate(triangle[:max_rows]):
        if len(row) < 2:
            break
        counter = Counter(row)
        total = len(row)
        max_val = max(row)
        mean_val = sum(row) / total

        # Fraction of 0s and 1s
        frac_01 = (counter.get(0, 0) + counter.get(1, 0)) / total

        if i < 15 or i % 10 == 0:
            print(
                f"\nRow {i:3d} (len={total:5d}): max={max_val:4d}, "
                f"mean={mean_val:.3f}, frac(0,1)={frac_01:.4f}"
            )
            # Show top-5 values
            top5 = counter.most_common(5)
            print(f"         Top values: {dict(top5)}")


def analyze_value_convergence(triangle):
    """
    KEY INSIGHT HUNT: Do the rows converge to a specific distribution?
    If rows quickly become dominated by 0s and 1s, that would explain
    why the first element stays 1.
    """
    print("\n" + "=" * 60)
    print("CONVERGENCE ANALYSIS: Do rows become dominated by {0,1}?")
    print("=" * 60)

    rows_data = []
    for i, row in enumerate(triangle):
        if len(row) < 2:
            break
        total = len(row)
        counter = Counter(row)
        frac_0 = counter.get(0, 0) / total
        frac_1 = counter.get(1, 0) / total
        frac_2 = counter.get(2, 0) / total
        max_val = max(row)
        rows_data.append(
            {
                "row": i,
                "frac_0": frac_0,
                "frac_1": frac_1,
                "frac_2": frac_2,
                "max": max_val,
                "len": total,
            }
        )

    print(
        f"\n{'Row':>5} {'%zeros':>8} {'%ones':>8} {'%twos':>8} {'%(0+1)':>8} {'max':>6}"
    )
    print("-" * 50)
    for d in rows_data[:30]:
        pct_01 = d["frac_0"] + d["frac_1"]
        print(
            f"{d['row']:5d} {d['frac_0']:8.4f} {d['frac_1']:8.4f} "
            f"{d['frac_2']:8.4f} {pct_01:8.4f} {d['max']:6d}"
        )

    return rows_data


# ============================================================
# PART 3: The "Almost-Binary" Hypothesis
# ============================================================


def test_binary_absorption(n_primes=500):
    """
    HYPOTHESIS: After a few rows, the sequence becomes "almost binary"
    (mostly 0s and 1s). When you take |a-b| of two values that are
    each 0 or 1, you get 0 or 1. So the sequence is self-reinforcing.

    The question becomes: WHY does the initial prime sequence lead to
    this almost-binary state so quickly?
    """
    print("\n" + "=" * 60)
    print("BINARY ABSORPTION HYPOTHESIS")
    print("=" * 60)
    print("If rows become mostly 0s and 1s, the process is self-reinforcing:")
    print("|0-0|=0, |0-1|=1, |1-0|=1, |1-1|=0")
    print("So a row of {0,1} maps to another row of {0,1}!")
    print()

    primes = [int(prime(i)) for i in range(1, n_primes + 1)]
    current = primes

    for i in range(min(n_primes - 1, 40)):
        next_row = [abs(current[j + 1] - current[j]) for j in range(len(current) - 1)]

        if len(next_row) < 2:
            break

        # Count how many elements are in {0, 1}
        binary_count = sum(1 for x in next_row if x <= 1)
        total = len(next_row)
        pct = binary_count / total * 100

        # Count "disruptions" (values > 1)
        disruptions = [
            (j, next_row[j]) for j in range(len(next_row)) if next_row[j] > 1
        ]
        n_disrupt = len(disruptions)

        if i < 15 or i % 5 == 0:
            print(
                f"Row {i + 1:3d}: {pct:6.2f}% binary ({binary_count}/{total}), "
                f"{n_disrupt} disruptions, max_disrupt={max(next_row) if next_row else 0}"
            )

        current = next_row


# ============================================================
# PART 4: What's Special About Primes vs Random Sequences?
# ============================================================


def compare_with_random(n_vals=200, n_trials=20):
    """
    CRITICAL TEST: Does Gilbreath's property hold for random odd-number
    sequences with similar gaps? If it does, it's about the STRUCTURE
    of the differences, not about primality per se. If it doesn't,
    then something specific about primes is essential.
    """
    print("\n" + "=" * 60)
    print("PRIMES vs RANDOM: Is primality essential?")
    print("=" * 60)

    # Test 1: Actual primes
    primes = [int(prime(i)) for i in range(1, n_vals + 1)]
    prime_result = check_gilbreath(primes)
    print(f"\nActual primes ({n_vals}): First column all 1s? {prime_result}")

    # Test 2: Random increasing odd integers with similar gaps
    import random

    random.seed(42)

    fail_count = 0
    for trial in range(n_trials):
        # Generate random increasing sequence starting with 2
        # with gaps similar to primes
        seq = [2]
        for _ in range(n_vals - 1):
            gap = random.choice([2, 4, 6, 2, 2, 4, 6, 8, 4, 2])  # prime-like gaps
            seq.append(seq[-1] + gap)

        result = check_gilbreath(seq)
        if not result:
            fail_count += 1

    print(
        f"Random prime-like sequences ({n_trials} trials): {fail_count}/{n_trials} FAILED"
    )

    # Test 3: Sequence starting with 2,3 then random odd numbers
    fail_count_2 = 0
    for trial in range(n_trials):
        seq = [2, 3]
        for _ in range(n_vals - 2):
            gap = random.choice([2, 4, 6, 8, 10, 12])
            seq.append(seq[-1] + gap)

        result = check_gilbreath(seq)
        if not result:
            fail_count_2 += 1

    print(
        f"Sequences starting 2,3 + random even gaps ({n_trials} trials): "
        f"{fail_count_2}/{n_trials} FAILED"
    )

    # Test 4: What if we start with 2,3 and force gap=2 (like twin primes everywhere)?
    twin_seq = [2] + [3 + 2 * i for i in range(n_vals - 1)]
    twin_result = check_gilbreath(twin_seq)
    print(
        f"All-twin-primes sequence (2,3,5,7,9,...): First column all 1s? {twin_result}"
    )


def check_gilbreath(seq):
    """Check if a sequence satisfies Gilbreath's property."""
    current = list(seq)
    for i in range(len(seq) - 1):
        next_row = [abs(current[j + 1] - current[j]) for j in range(len(current) - 1)]
        if len(next_row) == 0:
            break
        if next_row[0] != 1:
            return False
        current = next_row
    return True


# ============================================================
# PART 5: Disruption Propagation Analysis
# ============================================================


def track_disruptions(n_primes=300):
    """
    Track how 'large' values (>1) appear and propagate through rows.
    KEY QUESTION: Do disruptions (values > 1) always "heal" before
    reaching the first column?
    """
    print("\n" + "=" * 60)
    print("DISRUPTION PROPAGATION TRACKING")
    print("=" * 60)
    print("Tracking where values > 1 appear and how they move.\n")

    primes = [int(prime(i)) for i in range(1, n_primes + 1)]
    current = primes

    for i in range(min(n_primes - 1, 25)):
        next_row = [abs(current[j + 1] - current[j]) for j in range(len(current) - 1)]
        if len(next_row) < 2:
            break

        # Find positions of disruptions (values > 1)
        disruptions = [
            (j, next_row[j]) for j in range(min(len(next_row), 50)) if next_row[j] > 1
        ]

        # How close is the nearest disruption to position 0?
        if disruptions:
            nearest = min(j for j, v in disruptions)
            print(
                f"Row {i + 1:3d}: nearest disruption at pos {nearest:3d} "
                f"(value={next_row[nearest]}), total disruptions in first 50: {len(disruptions)}"
            )
        else:
            print(f"Row {i + 1:3d}: NO disruptions in first 50 elements! (pure binary)")

        current = next_row


# ============================================================
# MAIN: Run all analyses
# ============================================================

if __name__ == "__main__":
    print("🔬 GILBREATH'S CONJECTURE EXPLORER")
    print("=" * 60)
    print("Conjecture: In the prime difference triangle,")
    print("the first element of every row (after row 0) is always 1.")
    print("=" * 60)

    # Generate the triangle
    N = 500
    print(f"\nGenerating Gilbreath triangle for first {N} primes...")
    t0 = time.time()
    triangle = generate_gilbreath_triangle(N)
    print(f"Done in {time.time() - t0:.2f}s")

    # Verify conjecture
    print("\n--- VERIFICATION ---")
    valid, firsts = verify_conjecture(N)
    print(f"Conjecture holds for {N} primes: {valid}")
    print(f"First elements (row 0-9): {firsts[:10]}")

    # Print the triangle
    print("\n--- THE TRIANGLE (first 15 rows, 15 cols) ---")
    print_triangle(triangle, max_rows=15, max_cols=15)

    # Statistical analysis
    analyze_row_distributions(triangle, max_rows=30)

    # Convergence analysis
    rows_data = analyze_value_convergence(triangle)

    # Binary absorption
    test_binary_absorption(N)

    # Disruption tracking
    track_disruptions(N)

    # Compare with random
    compare_with_random(n_vals=200, n_trials=50)

    print("\n" + "=" * 60)
    print("EXPLORATION COMPLETE - See findings above!")
    print("=" * 60)
