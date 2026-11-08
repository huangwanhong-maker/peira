# Weighting

Status: **draft** (2026-11-08). Module `peira::weighting`.

## 1. Concepts

Weighting by the inverse of the probability of being seen: one estimator for samples and for arms. The theory, with proofs, is in [docs/theory.md](../docs/theory.md).

## 2. Content

- `horvitz_thompson_total`, `hajek_mean`, `ht_variance`, `syg_variance`.
- `ipw` (normalised or not), `aipw` (doubly robust, with its influence function), `trim`.
- `fit_nuisances`: the propensity by logistic regression and the arms' outcomes by least squares, through tyche.

## 3. Checks

Every sample of an unequal-probability design: the total and both variance estimators unbiased exactly; every assignment of a Bernoulli design: IPW unbiased exactly; double robustness by simulation.

## 4. Dependencies

peira `design`, tyche `probability` and `regression`, calx `linalg`.

## 5. References

Horvitz and Thompson (1952); Hájek (1971); Sen (1953), Yates and Grundy (1953); Särndal, Swensson and Wretman (1992); Robins, Rotnitzky and Zhao (1994); Hahn (1998).
