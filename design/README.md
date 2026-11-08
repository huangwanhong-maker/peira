# Designs

Status: **draft** (2026-11-08). Module `peira::design`.

## 1. Concepts

The assignment mechanisms by which units come to be treated, and the samples they need. The theory, with proofs, is in [docs/theory.md](../docs/theory.md).

## 2. Content

- `design.hpp`: Bernoulli, complete, blocked and clustered assignment, drawn from calx's streams; `each_complete` and `each_bernoulli` visit every assignment for exact randomisation distributions; `binomial_coefficient`.
- `power.hpp`: the power of a difference in means or proportions, and the sample per arm it needs, by the normal approximation; Kish's design effect.

## 3. Checks

Uniformity of complete randomisation (χ²), exact counts in blocks and clusters, Bernoulli enumeration; the 63 per arm for δ = σ/2 and 97 for proportions 0.6 against 0.4, and power against simulated z-tests.

## 4. Dependencies

calx `random`, tyche `probability`.

## 5. References

Imbens and Rubin (2015), *Causal Inference for Statistics, Social, and Biomedical Sciences*; Fisher (1935); Cohen (1988); Fleiss, Levin and Paik (2003); Kish (1965).
