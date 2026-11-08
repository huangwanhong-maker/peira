# Designs of Samples

Status: **draft** (2026-11-08). Module `peira::surveys`.

## 1. Concepts

Estimates from samples drawn by strata and clusters, every one a weighted total or a smooth function of totals, with the ultimate-cluster variance. The theory, with proofs, is in [docs/theory.md](../docs/theory.md).

## 2. Content

- `Design`: weights, strata, primary sampling units and sampling fractions.
- `total`, `ratio`, `mean`: estimates with standard errors on PSUs − strata degrees of freedom.
- `stratified_srswor`: the design of a stratified simple random sample.
- `rake`: iterative proportional fitting to known margins; post-stratification on one.

## 3. Checks

Every sample of stratified SRSWOR: the total unbiased, its variance the textbook formula, the ultimate-cluster estimator unbiased; every ordered sample of PSUs with replacement: Hansen–Hurwitz unbiased with an unbiased variance estimator; raking meets its margins.

## 4. Dependencies

peira `weighting`, tyche `probability`.

## 5. References

Wolter (2007), *Introduction to Variance Estimation*; Särndal, Swensson and Wretman (1992); Hansen and Hurwitz (1943); Woodruff (1971); Deming and Stephan (1940).
