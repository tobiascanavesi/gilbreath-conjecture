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
├── paper.tex                  # LaTeX source of the paper
├── requirements.txt           # Python dependencies
├── code/
│   ├── formal_argument.py     # Main formal argument with all propositions verified
│   ├── gilbreath_computation.py  # Comprehensive computations (bias, autocorrelation, triangle, |J_k|)
│   ├── hardy_littlewood_bias.py  # Hardy–Littlewood singular series and bias derivation
│   ├── gilbreath_explorer.py     # Exploratory analysis of the Gilbreath triangle
│   ├── mod4_analysis.py          # Mod 4 independence and correlation analysis
│   └── verify_results.py         # Quick verification of generated data files
├── data/
│   ├── bias_decay_data.csv       # Bias ε(x) at 8 scales (10²–10⁹)
│   ├── autocorrelation_data.csv  # Autocorrelation at lags 1–20
│   ├── xor_convergence_data.csv  # Theoretical XOR convergence curves
│   └── jk_size_data.csv          # |J_k| values for k = 1..1000
└── figures/
    ├── fig_heartbeat.pdf         # Gilbreath triangle heartbeat mechanism
    ├── fig_bias_decay.pdf        # Bias decay ε(x) vs log(x)
    ├── fig_xor_convergence.pdf   # XOR bias absorption convergence
    └── fig_autocorrelation.pdf   # Autocorrelation of b(n) sequence
```

## Quick start

```bash
pip install -r requirements.txt

# Run the main formal argument (verifies all propositions)
python code/formal_argument.py

# Run comprehensive computations and generate data files
python code/gilbreath_computation.py

# Explore the Gilbreath triangle interactively
python code/gilbreath_explorer.py

# Derive the Hardy–Littlewood bias
python code/hardy_littlewood_bias.py

# Analyze mod 4 independence
python code/mod4_analysis.py
```

## Key results

- **Bias decay**: ε(x) ≈ 0.97/log(x), confirmed at 8 scales from 10² to 10⁹
- **Weak mixing**: lag-1 autocorrelation ≈ −0.04; all lags ≥ 2 indistinguishable from zero
- **XOR absorption**: residual bias after m XORs is (1/2)(2ε)^m → 0 exponentially
- **Gilbreath verification**: d_k(1) = 1 verified for all k up to 5,000 rows

## Citation

If you use this code, please cite:

```bibtex
@article{canavesi2026gilbreath,
  title={Toward Gilbreath's Conjecture via XOR Bias Absorption of Prime Gap Parities},
  author={Canavesi, Tobias},
  year={2026}
}
```

## License

MIT License
