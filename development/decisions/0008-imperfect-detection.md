# 0008. Imperfect Detection: Occupancy, Repeated Counts, Capture–Recapture and Distance Sampling

- **Status:** Accepted, provisionally (2026-11-12).
- **Decided by:** the infrastructure's ADR 0101, which founds the ecology domain for urban forests and the life around them; the module's content, Claude.

## 1. Context

Surveys of organisms never see everything. A species present at a site may go undetected on a visit, an animal near a transect may not be seen, and a population is caught only in part. Ecology's first questions are about such surveys: where is it, how many are there, how dense are they. Each answer needs the probability of detection, which no design sets and which must be estimated from the data.

This is peira's structure (theory note, §1). The full data are the true occupancy $z_i$, the true abundance $N_i$, the population size $N$, or the individuals' positions. Detection coarsens them with probability $p$. Repeated visits, marks and distances are designs that identify the mechanism. The methods name no organism, and other sciences need them too: the survey of any object an observer can miss. By the infrastructure's rule they go into the foundation, and the ecology domain uses them.

## 2. Decision

A module `detection`, in the theory note's §10:
- `latent.hpp`: one likelihood for latent counts seen through binomial detection.
  - The mixing laws are Bernoulli, Poisson and negative binomial.
  - The gradient is exact, and so is the Hessian, by Louis's identity.
  - The fit is calx's BFGS polished by Newton's steps.
  - Occupancy is its Bernoulli case.
- `occupancy.hpp`: MacKenzie et al. (2002).
  - Logit-linear ψ and p, with standard errors and AIC.
  - The posterior of each site's occupancy.
  - The naive and finite-sample occupancy.
- `n_mixture.hpp`: Royle (2004), Poisson or negative binomial. The bound K is stated by the caller, and the fit reports two effects of it: the prior mass beyond K, and the change in the log-likelihood when K is doubled.
- `capture_recapture.hpp`:
  - Lincoln–Petersen, and Chapman with Seber's variance and Chao's interval;
  - Schnabel;
  - models M0 and Mt, by maximum likelihood over the integers N ≥ M, with profile likelihood intervals.
- `distance.hpp`: lines and points.
  - The keys are half-normal, hazard rate and uniform, with optional cosine adjustments.
  - The effective half-width or area has its delta-method cv.
  - The encounter rate takes Fewster et al.'s R2 variance.
  - Density has Satterthwaite's degrees of freedom and a log-normal interval.
  - AIC chooses among detection functions.
- `simulate.hpp`: surveys drawn from known truths.

The theory note's §10 states the structure the four designs share. Each estimate is a Horvitz–Thompson sum over what was seen, with an estimated probability of being seen at least once. Each design replicates the mechanism so that this probability is identified. Distance sampling is `weighting`'s estimator, and a check holds it so.

## 3. Alternatives Considered

- **In the ecology domain.** Rejected: the methods name no organism, and epidemiology and any other survey need them.
- **In tyche, as latent-variable likelihoods.** Rejected: their object is the mechanism of observation and what it identifies, which is peira's boundary (decision 0002). tyche keeps inference within a given model.
- **Separate likelihoods for occupancy and the N-mixture model.** Rejected: occupancy is the N-mixture model with a latent count of 0 or 1. One sum, with Louis's identity, gives both their exact information.
- **Hessians by differences for the latent models.** Rejected where an exact one is cheap. For distance sampling, whose models have one to four parameters and integrals by quadrature, the Hessian is by differences.
- **Deferred:**
  - heterogeneity among individuals (models $M_h$, $M_b$), whose $N$ is identified only within a family (Link 2003);
  - open populations;
  - occupancy across seasons, and spatial capture–recapture;
  - covariates in the scale of a detection function, clusters, and grouped distances.

## 4. Verification

Four check programs, 124 checks. Simulations are the backbone, beside closed forms and enumeration:
- **Likelihoods.** Each equals a brute-force or written-out sum on small data, and the gradients and Hessians agree with differences.
- **Identification.** With one visit, occupancy depends on ψp alone and the N-mixture model on λp alone (thinning). A visit made by no site changes nothing.
- **Closed forms:**
  - the constant occupancy model's characterisation by the truncated binomial;
  - Chapman's estimator exactly unbiased when n1 + n2 ≥ N, with its exact bias otherwise (the theory note's Theorem 10.5, proved there);
  - M0 and Mt equal to exhaustive searches over N;
  - the half-normal's σ̂², μ and cv without truncation;
  - R2 by hand.
- **Horvitz–Thompson.** D̂ × the covered area equals `weighting`'s total.
- **Recovery.** Surveys simulated from known truths, for occupancy, the N-mixture model, M0, Mt, and line and point transects, recover the truth:
  - bias below a quarter (0.3 for the N-mixture model and points) of the estimates' spread;
  - coverage within 3.5 Monte Carlo standard errors of 95%;
  - standard errors within 15–20% of the spread.

No reference package was run: R's `unmarked` and `Distance` are not installed.

## 5. Consequences

- The infrastructure's ecology domain composes these designs for urban forests and the life around them.
- peira now depends on calx's `optimisation`, `quadrature` and `special_functions`, through `detection` alone.
- peira's checks rise from 168 to 292.
