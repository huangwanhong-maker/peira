# 0003. The First Modules

- **Status:** Accepted, provisionally
- **Date:** 2026-11-08
- **Decided by:** Claude, after decision 0002 and the stages of the infrastructure's ADR 0082

## 1. Context

The first stage should show the structure of decision 0002 in its four forms: designs that assign, designs that sample, censoring and nonresponse. Each module needs a truth to be checked against.

## 2. Decision

Six modules, in this order:

| Module | Content | Checked against |
|---|---|---|
| `design` | Assignment mechanisms (Bernoulli, complete, blocked, clustered), their enumeration, power and sample size | uniformity of draws; closed forms of power |
| `randomisation` | Neyman's difference in means and its conservative variance, the stratified estimator, Fisher's randomisation test | every assignment of small populations: unbiasedness, Neyman's variance identity, the level of the test; Fisher's lady tasting tea |
| `weighting` | Horvitz–Thompson, Hájek, their variance estimators; inverse-probability weighting and the doubly robust estimator; the nuisance models through tyche | every sample of an unequal-probability design; every assignment of a Bernoulli design; simulated double robustness |
| `surveys` | Totals, means and ratios under designs of strata and clusters with the ultimate-cluster variance; raking | every sample of stratified SRSWOR and of PSUs drawn with replacement |
| `survival` | Kaplan–Meier, Nelson–Aalen, log-rank, Cox (Breslow, Efron), Aalen–Johansen | statsmodels 0.15.0 and lifelines 0.30.3 on Rossi's recidivism study; exact identities |
| `missing` | Rubin's rules with Barnard–Rubin degrees of freedom; the approximate Bayesian bootstrap | closed forms; simulated coverage under MAR |

Each is a header-only target `peira::<module>` with a check program `check_peira_<module>`.

## 3. Alternatives Considered

- **Quasi-experiments first.** They wait for the second stage, so that they move from nomos with their fixed effects going to tyche first, rather than being written twice.
- **Cohorts (person-time, standardisation, age–period–cohort) in the first stage.** They follow in the third, with the Lexis diagram's arithmetic, once survival's risk sets are in place.

## 4. Consequences

- 89 checks in the first stage.
- The next stages, in the infrastructure's ADR 0082:
  - the quasi-experiments and causal graphs;
  - then cohorts, meta-analysis, sequential testing, sensitivity and pre-registered plans.
