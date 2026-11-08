# 0004. Quasi-Experiments and Causal Graphs: The Second Stage

- **Status:** Accepted, provisionally
- **Date:** 2026-11-08
- **Decided by:**
  - **The plan.** The infrastructure's ADR 0082, stage 2: differences in differences and instruments move from nomos, panel fixed effects to tyche.
  - **The go-ahead.** The user ("all proceed", 2026-11-08).
  - **The modules.** Claude.

## 1. Context

Stage 1 (decision 0003) held the designs whose mechanism is known. Stage 2 holds identification by assumptions on the outcome (parallel trends, exclusion, continuity) and by a causal graph (theory note, §7–8). nomos held differences in differences and instruments because the economics domain needed them first; they serve every science that observes.

## 2. Decision

Two modules.

### 2.1 `quasi`

- **Moved from nomos, with their checks.**
  - `difference_in_differences.hpp`: the 2 × 2 estimator and the two-way regression;
  - `instrumental.hpp`: 2SLS, the first stage, Sargan and Hansen.

  Their panels and robust covariances are now tyche's (`tyche::panel`, `tyche::regression`'s sandwich; tyche's decision 0014).
- **New.**
  - `staggered.hpp`: Callaway and Sant'Anna's group-time effects against never or not-yet treated units, aggregated simply and by exposure, with a unit bootstrap; and the static and dynamic two-way regressions beside them.
  - `instrumental.hpp`'s Wald estimator of the local average effect.
  - `discontinuity.hpp`: the sharp and fuzzy designs by local linear regression, three kernels, a cross-validated bandwidth at the boundary.
  - `synthetic.hpp`: simplex weights by accelerated projected gradient, with permutation inference in space.

### 2.2 `graphs`

`causal.hpp`: directed acyclic graphs, d-separation by reachability, the back-door and front-door criteria, and minimal adjustment sets by enumeration. It is written here rather than on morphe's graphs, as a dozen lines of reachability need no more.

## 3. Checked Against

- **Moved checks:** nomos's checks of the 2 × 2 design and the instruments, unchanged.
- **Staggered adoption:** noise-free panels whose effects differ by group and grow with exposure. Every group-time effect is recovered exactly, the placebos are zero, the two-way regression's departure (1.94 against 2.43) is shown, and the two agree under one constant effect.
- **The Wald estimator:** it is 2SLS exactly, and recovers the compliers' effect.
- **Regression discontinuity:** lines are exact at any bandwidth; the bias falls fourfold when the bandwidth halves; the fuzzy design is the ratio of the jumps.
- **Synthetic control:** exact convex combinations give back their weights; the optimality conditions on the simplex hold.
- **Causal graphs:** d-separation against vanishing partial correlations in linear Gaussian models, 4,032 cases of 4,032.

## 4. Consequences

- 48 checks more (137 in all).
- nomos keeps marginal effects, the multinomial logit and Poisson regression with absorbed effects (its decision 0007).
- The API's toolbox reads panels from tyche and instruments and differences in differences from peira.
- Next, stage 3: `cohorts` (person-time, rates, standardisation, the Lexis diagram, age–period–cohort), on which the demography domain stands.
