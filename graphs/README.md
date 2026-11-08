# Causal Graphs

Status: **draft** (2026-11-08). Module `peira::graphs`.

## 1. Concepts

A directed acyclic graph says which variables cause which. It implies conditional independences (d-separation), and tells which observed sets identify an effect by adjustment: the back-door and front-door criteria. The theory is in [docs/theory.md](../docs/theory.md), §8.

## 2. Content

- `causal.hpp`:
  - `Dag`, `acyclic`, `descendants`;
  - `d_separated` (reachability of active trails);
  - `without_outgoing`, `backdoor`, `frontdoor`;
  - `adjustment_sets`, the minimal back-door sets among the observed.

## 3. Checks

`check_peira_graphs` (11):
- chains, forks and colliders;
- M-bias;
- the front door with a hidden confounder;
- d-separation against the vanishing of partial correlations in random linear Gaussian models (4,032 cases);
- acyclicity.

## 4. Dependencies

None in its header. In the checks, calx `linalg` (the covariance of a linear model) and `random`.

## 5. References

Pearl (1988, 2009); Verma and Pearl (1988); Geiger, Verma and Pearl (1990); Shachter (1998); Koller and Friedman (2009); Spirtes, Glymour and Scheines (2000).
