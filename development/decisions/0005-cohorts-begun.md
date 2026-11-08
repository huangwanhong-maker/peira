# 0005. Cohorts Begun: the Lexis Diagram, Rates and Standardisation

- **Status:** Accepted, provisionally (2026-11-08).
- **Decided by:** the infrastructure's ADR 0085 (demography founded on peira, at the user's word); the module's content, Claude.

## 1. Decision

The third stage begins with `cohorts`, what the demography domain needs first:
- `cohorts.hpp`: lives split on the Lexis diagram into person-time and events by age and period, exactly, with what falls outside reported;
- rates with Garwood's exact Poisson intervals;
- direct standardisation with its Poisson variance, and indirect standardisation (the SMR) with Garwood's interval.

Age–period–cohort analysis, with its non-identification shown exactly, follows in this module.

## 2. Verification

`check_peira_cohorts` (9):
- one life split by hand;
- conservation of person-time and events over 400 random lives;
- agreement with a brute-force integration;
- Garwood's bounds against the Poisson's tails (to 10⁻¹⁶);
- standardisation's identities.

## 3. Consequences

The demography domain reads its death rates from the Lexis diagram, and its standardised rates through `direct`.
