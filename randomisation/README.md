# Randomisation Inference

Status: **draft** (2026-11-08). Module `peira::randomisation`.

## 1. Concepts

Inference for randomised experiments from the design alone: potential outcomes are fixed, and the assignment is the only randomness. The theory, with proofs, is in [docs/theory.md](../docs/theory.md).

## 2. Content

- `difference_in_means`: Neyman's estimator and its conservative variance, with a normal interval.
- `stratified`: blocks' estimates weighted by their shares.
- `fisher`: the randomisation test of a sharp null of constant effect τ₀, enumerated when it is small enough, else drawn; one- or two-sided.

## 3. Checks

Every assignment of small populations: E[τ̂] = τ and Neyman's variance identity exactly, conservative with heterogeneous effects and exact with a constant one; the lady tasting tea (p = 1/70); the level of Fisher's test over all 252 assignments; coverage over simulated experiments.

## 4. Dependencies

peira `design`, tyche `probability` and `testing`.

## 5. References

Neyman (1923); Fisher (1935); Imbens and Rubin (2015), ch. 5–6, 9; Phipson and Smyth (2010).
