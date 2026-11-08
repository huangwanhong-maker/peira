# Cohorts

Status: **draft** (2026-11-08). Module `peira::cohorts`.

## 1. Concepts

A cohort is seen in time: each life a segment on the Lexis diagram of calendar time and age, cut by its grid into person-time and events. Rates are events over person-time, and rates of populations of different structure are compared by standardising them. The general methods of demography and epidemiology begin here; the theory is in [docs/theory.md](../docs/theory.md), §9, and the module's header.

## 2. Content

- `cohorts.hpp`:
  - `Life`, `Grid`, `lexis`: person-time and events by age and period band, split exactly;
  - `rate`, with Garwood's exact interval;
  - `direct` and `indirect` standardisation.

## 3. Checks

`check_peira_cohorts` (9):
- a life split by hand;
- the conservation of person-time and events;
- agreement with brute-force integration;
- Garwood's bounds against the Poisson's tails;
- standardisation's identities.

## 4. Dependencies

tyche `probability` (the χ² quantile, the normal).

## 5. References

Lexis (1875); Keiding (1990); Garwood (1936); Breslow and Day (1987, ch. 2); Fay and Feuer (1997).
