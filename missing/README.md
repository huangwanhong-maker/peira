# Missing Data

Status: **draft** (2026-11-08). Module `peira::missing`.

## 1. Concepts

Inference when values are missing at random, by several imputations combined by Rubin's rules. The theory, with proofs, is in [docs/theory.md](../docs/theory.md).

## 2. Content

- `pool`: Rubin's rules, with Barnard and Rubin's degrees of freedom and the fraction of information missing.
- `approximate_bayesian_bootstrap`: a proper imputation within classes.

## 3. Checks

Rubin's rules in closed form; under MAR, the complete-case mean biased and multiple imputation covering at its nominal rate in simulation.

## 4. Dependencies

peira `design`, tyche `probability`.

## 5. References

Rubin (1976, 1987); Rubin and Schenker (1986); Barnard and Rubin (1999); van Buuren (2018), *Flexible Imputation of Missing Data*.
