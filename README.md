# Toward Gilbreath's Conjecture via XOR Bias Absorption

Companion code and data for the paper:

> **Toward Gilbreath's Conjecture via XOR Bias Absorption of Prime Gap Parities**
> Tobias Canavesi, March 2026

## Overview

Gilbreath's conjecture (1958) states that the first element of every row in the iterated absolute-difference triangle of primes is always 1. This repository contains the computational framework that supports a conditional proof via three mechanisms:

1. **Mod 4 reduction** (Odlyzko): the Gilbreath triangle mod 4 follows XOR dynamics on binary sequences derived from prime gaps.
2. **Hardy–Littlewood bias**: these binary sequences carry a computable bias ε(x) ≈ c/log(x) from the singular series.
3. **XOR bias absorption** (new): XOR sums of weakly-dependent biased bits converge exponentially to fairness.

## Repository structure

```
├── paper.tex                     # LaTeX source of the paper
├── requirements.txt              # Python dependencies
├── code/
│   ├── formal_argument.py        # Main formal argument with all propositions verified
│   ├── gilbreath_computation.py  # Core computations (bias, autocorrelation, triangle, |J_k|)
│   ├── gilbreath_extended.py     # Extended computations up to 10^14 (bias), 10^10 (autocorr)
│   ├── hardy_littlewood_bias.py  # Hardy–Littlewood singular series and bias derivation
│   ├── gilbreath_explorer.py     # Exploratory analysis of the Gilbreath triangle
│   ├── mod4_analysis.py          # Mod 4 independence and correlation analysis
│   └── verify_results.py        # Quick verification of generated data files
├── data/
│   ├── bias_decay_data.csv             # Bias ε(x) at 8 scales (10²–10⁹)
│   ├── bias_decay_extended.csv         # Bias ε(x) at 13 scales (10²–10¹⁴), 200k primes
│   ├── autocorrelation_data.csv        # Autocorrelation at lags 1–20 (4 scales)
│   ├── autocorrelation_extended.csv    # Autocorrelation at lags 1–50 (8 scales, up to 1M primes)
│   ├── absorption_threshold_extended.csv  # k* vs G_N at 10 scales (N up to 500k)
│   ├── absorption_convergence.csv      # % entries in {0,2} at each row (200k primes)
│   ├── xor_convergence_data.csv        # Theoretical XOR convergence curves
│   └── jk_size_data.csv               # |J_k| values for k = 1..1000
└── figures/
    ├── fig_heartbeat.pdf         # Gilbreath triangle heartbeat mechanism
    ├── fig_bias_decay.pdf        # Bias decay ε(x) vs log(x)
    ├── fig_xor_convergence.pdf   # XOR bias absorption convergence
    └── fig_autocorrelation.pdf   # Autocorrelation of b(n) sequence
```

## Reproducing the results

```bash
pip install -r requirements.txt

# Run the main formal argument (verifies all propositions)
python code/formal_argument.py

# Run core computations and generate data files (Table 1, Figures 1–4)
python code/gilbreath_computation.py

# Run extended computations — up to 10^14, 1M primes per scale (~3 min)
python code/gilbreath_extended.py

# Explore the Gilbreath triangle interactively
python code/gilbreath_explorer.py

# Derive the Hardy–Littlewood bias from the singular series
python code/hardy_littlewood_bias.py

# Analyze mod 4 independence structure
python code/mod4_analysis.py
```

## Key results

### Bias decay (ε(x) ≈ c/log x)
Measured at 13 scales from 10² to 10¹⁴ with up to 200,000 primes per scale.
The product ε·log(x) = **0.984 ± 0.034**, confirming ε(x) ~ 1/log(x).
Power-law fit exponent α ≈ 0.94, converging toward the conjectured α = 1.

### 2-dependence of prime gap parities
Measured at 8 scales from 10³ to 10¹⁰ with up to 1,000,000 primes per scale.
- Lag-1 autocorrelation: −0.04 to −0.02 (significant, decaying with scale)
- Lag-2 autocorrelation: +0.006 to +0.002 (significant at smaller scales, vanishing)
- Lag ≥ 3: never systematically significant (consistent with chance at 5% level)

### Binary absorption threshold
k*/G_N = **0.825 ± 0.091** across 10 scales (N = 500 to 500,000), supporting Conjecture 7.2.
Full convergence-to-absorption curve provided: 10% at row 1, 88% at row 20, 99.2% at row 50, 100% at row 113 (for 200k primes).

### Gilbreath verification
d_k(1) = 1 verified for **all 2,000 rows** with 200,000 primes. Zero violations.
Reset frequency: 49.1%, max gap between resets: 12.

## Citation

```bibtex
@article{canavesi2026gilbreath,
  title={Toward Gilbreath's Conjecture via XOR Bias Absorption of Prime Gap Parities},
  author={Canavesi, Tobias},
  year={2026}
}
```

## License

MIT License
