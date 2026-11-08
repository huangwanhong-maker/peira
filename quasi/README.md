# Quasi-Experiments

Status: **draft** (2026-11-08). Module `peira::quasi`.

## 1. Concepts

Where the mechanism of assignment is unknown, an assumption on the potential outcomes identifies an effect:
- **parallel trends:** differences in differences, and staggered adoption by group and period;
- **exclusion and monotonicity:** instruments, and the local average effect;
- **continuity at a cut-off:** regression discontinuity;
- **a past reproduced by a weighted average of others:** synthetic control.

The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §7.

## 2. Content

- `difference_in_differences.hpp` (moved from nomos): the 2 × 2 estimator with clustered errors; the two-way fixed-effects regression.
- `staggered.hpp`: ATT(g, t) against never or not-yet treated units (Callaway and Sant'Anna 2021), their simple and exposure aggregations with a unit bootstrap; the static and dynamic two-way regressions for comparison.
- `instrumental.hpp` (moved from nomos): two-stage least squares with classical, HC0, HC1 and clustered covariances, the first stage's F and partial R², Sargan's and Hansen's tests; `wald`, the local average treatment effect with a binary instrument.
- `discontinuity.hpp`: sharp and fuzzy regression discontinuity by local linear regression (triangular, Epanechnikov, uniform kernels), robust standard errors, a bandwidth by leave-one-out cross-validation at the boundary.
- `synthetic.hpp`: synthetic control weights on the simplex (accelerated projected gradient, the exact projection), gaps, RMSPE before and after, placebo inference in space; `project_simplex`.

## 3. Checks

`check_peira_quasi` (26) and `check_peira_instrumental` (11, moved from nomos). The truths are exact:
- noise-free panels and designs;
- the identities of the Wald and fuzzy estimators with 2SLS;
- the boundary bias's rate as the bandwidth halves;
- simplex optimality;
- coverage within four standard errors with noise.

## 4. Dependencies

tyche `panel` (within estimators, absorbed effects), `regression` (least squares, the sandwich), `probability`; calx `linalg`, `random`.

## 5. References

Card and Krueger (1994); Callaway and Sant'Anna (2021); Goodman-Bacon (2021); Sun and Abraham (2021); Theil (1953); Imbens and Angrist (1994); Angrist, Imbens and Rubin (1996); Hahn, Todd and van der Klaauw (2001); Imbens and Lemieux (2008); Ludwig and Miller (2007); Abadie, Diamond and Hainmueller (2010); Duchi et al. (2008); Beck and Teboulle (2009).
