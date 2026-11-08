# 0002. Scope, Boundaries and Dependencies

- **Status:** Accepted
- **Date:** 2026-11-08
- **Decided by:**
  - **The request.** The user proposed "a quantitative (experiment) science foundation" for "classical and modern advanced methods for social science experiment, cohort analysis & visualization, hypothesis validation", observing that such methods serve "economics, clinical experiments, social science, psychology". The user agreed to the name peira, to a scope including surveys, survival and missing data, and to moving difference in differences and instrumental variables out of nomos.
  - **The criterion.** Claude, after the infrastructure's ADR 0082.

## 1. Context

The infrastructure placed its first methods of causal inference where the first science that needed them stood: difference in differences and instrumental variables in nomos, because the economics domain needed them. But a clinical trial, a psychology experiment, a survey of households and a cohort of patients use the same methods. Methods that name no object in nature belong in a foundation, and these did not belong to the social sciences.

What they share is one structure. In each case the full data are hidden in part by a mechanism, and the question is what can still be identified through it:
- in an experiment, assignment reveals one potential outcome of each unit;
- in a survey, sampling reveals some units;
- in a cohort, censoring reveals part of each life;
- with missing values, nonresponse reveals some values.

The Horvitz–Thompson estimator of a survey and the inverse-probability-weighted estimator of an experiment are the same estimator. Rubin's fundamental problem of causal inference is a missing-data problem. Kaplan–Meier is an inverse-probability-of-censoring estimator. The theory note (`docs/theory.md`) states this structure.

## 2. Decision

### 2.1 Scope

peira holds the designs and mechanisms by which units come to be treated and observed, and what is identified through them:
- potential outcomes and their estimands;
- randomised designs and their power;
- randomisation inference;
- weighting by the probability of being seen;
- designs of samples and their estimates;
- survival and censoring, and cohorts;
- missing data and multiple imputation;
- quasi-experiments: difference in differences, instrumental variables, regression discontinuity, synthetic control;
- adjustment sets on causal graphs;
- meta-analysis and sequential testing;
- sensitivity analysis and partial identification;
- pre-registered analysis plans.

**The criterion.** A method belongs here when its object is the mechanism: the law of assignment, of sampling, of censoring or of response, and what it lets one identify. A method that infers within a given model belongs to tyche. A measure or mechanism of social life belongs to nomos.

| Topic | Where it belongs |
|---|---|
| Distributions, likelihood, estimators and their asymptotics, regression and generalised linear models, panels with fixed effects, multiple testing, permutation machinery, time series | tyche |
| Inequality, index numbers, games, matching markets, voting, agents | nomos |
| Causal graphs as graphs: d-separation's combinatorics, shortest paths | morphe |
| Random streams, linear algebra | calx |
| A population, a disease, a trial of a particular drug, a city's households | the domains of the infrastructure |
| Drawing a Kaplan–Meier curve, a forest plot, a Lexis diagram | the infrastructure's `visualization/` |

### 2.2 Dependencies

peira depends on calx, morphe, physis and tyche, through their public headers, and on nothing else of the infrastructure. It does not depend on nomos. nomos may come to depend on peira, if one of its methods needs a design.

### 2.3 Invariants

1. **Pure.** No I/O, no global state; every draw from a calx stream given by the caller or seeded explicitly.
2. **Truth from the design.** Wherever the design can be enumerated, a check enumerates it, and the expectation is exact: unbiasedness, variance identities, the level of a test. Otherwise a reference library's values, made by running it (a script in `tests/fixtures/`), or a simulation whose Monte Carlo error is stated.
3. **Identification stated.** Each estimator's header says what it identifies and under which assumptions; assumptions that cannot be tested are said to be so.
4. **Cited sources.** Every method cites its source; equation numbers only where certain.

## 3. Alternatives Considered

- **A domain of "quantitative social science".** Rejected: a domain holds the objects of a science, and these are methods. The social sciences' domains (demography, sociology, political science, …) stand on peira.
- **Extending tyche.** tyche infers given a model; designs are about how the data came to be, whose theorems are of identification rather than estimation. Keeping them apart keeps both libraries' scopes crisp.
- **Leaving the methods in nomos.** Rejected, as the user observed: clinical and psychological experiments would then depend on a library of the social sciences.

## 4. Consequences

- The quasi-experiments move from nomos to peira, and panel fixed effects with Poisson regression with absorbed effects move from nomos to tyche, so that peira need not depend on nomos (ADR 0082's second stage).
- The infrastructure's map of mathematics places MSC 62D (sampling, missing data, causal inference from observational studies), 62K (design of experiments) and 62N (survival and censored data) here.
