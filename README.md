# peira

*πεῖρα, a trial; whence "empirical".* peira holds the designs and the mechanisms by which units come to be treated and observed, and what can be identified through them: experiments and their randomisation, samples and their designs, cohorts seen through censoring, values missing at random. It is a foundation library of [Noetarium, the Science and Humanity Infrastructure](https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure), and makes sense without it. Status: **draft** (2026-11-08).

## 1. Concepts

Every method here studies one structure ([docs/theory.md](docs/theory.md), §1). Full data are hidden in part by a mechanism, and the question is what the seen part still identifies:
- **Experiments:** assignment reveals one potential outcome of each unit.
- **Surveys:** sampling reveals some units.
- **Cohorts:** censoring reveals part of each life.
- **Missing data:** nonresponse reveals some values.

The Horvitz–Thompson estimator of a survey and the inverse-probability-weighted estimator of an experiment are one estimator. Fisher's exact test and Neyman's variance come from the design alone. Kaplan–Meier is a weighting by the probability of remaining uncensored. Rubin's rules carry the uncertainty of what nonresponse hid.

A method belongs here when its object is the mechanism and what it identifies. Inference within a given model is tyche's, and measures of social life are nomos's ([decision 0002](development/decisions/0002-scope-boundaries-and-dependencies.md)).

## 2. Modules

| Module | Holds | Status |
|---|---|---|
| [design](design/README.md) | Assignment mechanisms (Bernoulli, complete, blocked, clustered) and their enumeration; power and sample size | draft |
| [randomisation](randomisation/README.md) | Neyman's difference in means and its conservative variance; stratified estimates; Fisher's randomisation test | draft |
| [weighting](weighting/README.md) | Horvitz–Thompson, Hájek, their variance estimators; inverse-probability weighting; the doubly robust estimator | draft |
| [surveys](surveys/README.md) | Totals, means and ratios under strata and clusters, with the ultimate-cluster variance; raking and post-stratification | draft |
| [survival](survival/README.md) | Kaplan–Meier, Nelson–Aalen, the log-rank test, Cox's model (Breslow, Efron), Aalen–Johansen | draft |
| [missing](missing/README.md) | Rubin's rules with Barnard–Rubin degrees of freedom; the approximate Bayesian bootstrap | draft |
| [quasi](quasi/README.md) | Differences in differences (2 × 2 and two-way, moved from nomos), staggered adoption by group and period (Callaway and Sant'Anna) with event studies; instrumental variables (moved from nomos) and the Wald estimator of the local average effect; regression discontinuity, sharp and fuzzy, with a cross-validated bandwidth; synthetic control with inference in space | draft |
| [graphs](graphs/README.md) | Causal graphs: d-separation, the back-door and front-door criteria, minimal adjustment sets | draft |
| [cohorts](cohorts/README.md) | Person-time and events on the Lexis diagram, rates with exact intervals, direct and indirect standardisation (begun for demography); age–period–cohort and what it cannot identify, planned | draft |
| synthesis | Meta-analysis: fixed and random effects, heterogeneity, publication bias | planned |
| sequential | Group-sequential boundaries, alpha spending, always-valid inference | planned |
| sensitivity | Rosenbaum's bounds, E-values, Manski's bounds | planned |
| protocol | A pre-registered analysis plan as a declared object, and the deviations from it | planned |

## 3. Checks

`node development/tools/run-cpp.mjs foundations/peira foundations/peira` from the infrastructure's root builds and runs the 146 checks. Wherever a design can be enumerated, a check enumerates it: unbiasedness, variance identities and the level of a test then hold exactly. Survival is held to statsmodels 0.15.0 and lifelines 0.30.3 on Rossi's recidivism study (`survival/tests/fixtures/`, with the script that made them).

## 4. Dependencies

calx (random streams, linear algebra), morphe and physis (as the siblings' order requires) and tyche (distributions, regression, generalised linear models, panels with absorbed effects, permutation p-values). Header-only C++20.

## 5. Licence

MIT; see [LICENSE](LICENSE) and [CONTRIBUTORS.txt](CONTRIBUTORS.txt).
