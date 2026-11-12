# Detection

Status: **draft** (2026-11-12). Module `peira::detection`, [decision 0008](../development/decisions/0008-imperfect-detection.md).

## 1. Concepts

A survey never sees everything that is there. A site that holds something may yield no detection on a visit; an individual near a line may pass unseen; a population is caught only in part. The full data are what is there: whether a site is occupied, how many individuals a site or a population holds, where individuals lie. The mechanism is detection, which reveals each thing present with a probability $p < 1$ and never reveals what is absent. Unlike a sample drawn by design (`weighting`, `surveys`), the probability of being seen is not known. A single look cannot separate absent from missed, since nothing is seen with probability $(1 - \psi) + \psi(1 - p)$, a function of $\psi p$ alone.

Each design here replicates the mechanism so that its law can be identified:
- **Occupancy:** repeated visits while what a site holds stays fixed.
- **The N-mixture model:** repeated counts while what a site holds stays fixed.
- **Capture–recapture:** marks, which show who is caught again.
- **Distance sampling:** distances, with everything on the line or at the point assumed seen.

Every estimate is then a Horvitz–Thompson sum $\sum_{\text{seen}} 1/\hat\pi$ over what was seen, with $\hat\pi$ an *estimated* probability of being seen at least once:
- **occupancy:** $s\hat\psi = s_D / \hat p^*$, $p^* = 1 - (1 - p)^T$;
- **capture–recapture:** Darroch's equation $N = M/(1 - \prod_j (1 - \hat p_j))$;
- **distance sampling:** $\hat D \cdot 2wL = n/\hat P_a$, $P_a = \mu/w$ (this is `weighting`'s estimator, checked to $10^{-15}$).

Occupancy is also the N-mixture model whose latent count is 0 or 1, and one likelihood serves both. Each design's identifying assumption is stated in its header, and so is the part of it that the data cannot test. The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §10.

The module names no organism: its units are sites, visits, occasions, individuals, lines and points, and the ecology of the infrastructure, epidemiology or the survey of any object that can be missed composes them.

## 2. Content

- `common.hpp`:
  - `SiteVisits`, a table of sites by visits (counts, −1 for a visit not made);
  - `intercept`;
  - `FitOptions`;
  - `predict`: a probability or a mean at a row of covariates, with the delta method's standard error and an interval carried through the link.
- `latent.hpp`: latent counts seen through binomial detection.
  - The mixing law is `Mixing::bernoulli`, `poisson` or `negative_binomial`.
  - The likelihood sums over each site's count, in log space.
  - The gradient is exact, and so is Louis's (1982) Hessian.
  - The fit is calx's BFGS, then Newton's steps with the exact Hessian.
- `occupancy.hpp`: `occupancy` and `occupancy_log_likelihood` (MacKenzie et al. 2002, 2018).
  - Logit-linear ψ on site covariates and p on site-by-visit covariates.
  - Standard errors from the observed information, and AIC.
  - The posterior P(z_i = 1 | y_i).
  - The naive occupancy, and the finite-sample occupancy.
  - `predict_psi`, `predict_p`.
- `n_mixture.hpp`: `n_mixture` and `n_mixture_log_likelihood` (Royle 2004).
  - Poisson or negative-binomial abundance, log-linear λ and logit-linear p.
  - The bound K is stated by the caller, and its effect is reported: the prior mass beyond K, and the change in ℓ when K is doubled.
  - The posterior law of each N_i and its mean; the predicted and expected totals.
  - `predict_lambda`, `predict_p`.
- `capture_recapture.hpp`:
  - `lincoln_petersen`, and `chapman` with Seber's variance, both with Chao's (1987) log-normal interval;
  - `CaptureHistories`, and `summarise`;
  - `schnabel`, with the interval from Σ R as Poisson;
  - `closed_m0` and `closed_mt` (Otis et al. 1978; Darroch 1958): maximum likelihood over the integers N ≥ M, with profile likelihood intervals and AIC, from histories or from sufficient statistics;
  - `closed_log_likelihood`.
- `distance.hpp` (Buckland et al. 2001):
  - `DistanceSurvey`, for lines or points;
  - `DetectionModel`: a key, uniform, half-normal or hazard rate, with optional cosine adjustments;
  - `fit_detection`: the maximum likelihood estimate; the effective strip half-width μ or effective area ν with its cv by the delta method; f(0) and P_a;
  - `density`: the encounter rate with Fewster et al.'s (2009) R2 variance, D̂ with its cv, Satterthwaite's degrees of freedom and the log-normal interval, and the Horvitz–Thompson abundance in the covered area;
  - `distance_sampling`: several detection models compared by AIC;
  - `detection_probability`.
- `simulate.hpp`: surveys drawn from known truths, for planning and for checks.
  - `simulate_visits` (presence or counts);
  - `simulate_captures`;
  - `simulate_distances` (lines or points).
- `detection.hpp`: all of them.

Units: distances and lengths in one unit; density per that unit squared.

## 3. Checks

Four programs, 124 checks. Simulation is the backbone: every design is recovered from surveys drawn from a known truth. The other checks are closed forms, identities, and enumerations of likelihoods on small data.

- `check_peira_detection`, occupancy (28):
  - the likelihood equals the brute-force sum over all 2⁵ latent occupancies;
  - the gradient and Louis's Hessian agree with differences;
  - with one visit, ℓ depends on ψp alone, and a visit made by no site changes nothing;
  - the constant model's MLE is the truncated binomial's characterisation, solved by calx's Brent method;
  - the posterior of a site with no detection is Bayes's rule;
  - 500 simulated surveys of 300 sites: bias, Wald coverage and standard errors, the finite-sample occupancy unbiased for the realised, and the naive estimate short of it; AIC.
- `check_peira_n_mixture` (35):
  - the likelihood equals the sum written out with tyche's Poisson and binomial probabilities, and the negative binomial's;
  - the gradient and Hessian agree with differences;
  - thinning: one visit gives Poisson(λp), and the posterior of N − y is Poisson(λ(1 − p));
  - the negative binomial tends to the Poisson;
  - the bound K: negligible when reported so, too small when it sits at the largest count;
  - 200 simulated surveys: bias, coverage and standard errors, and the predicted total unbiased;
  - the negative binomial's r recovered on 1,000 sites.
- `check_peira_capture_recapture` (29):
  - Lincoln–Petersen and Chapman with Seber's variance, by hand;
  - Chapman's estimator exactly unbiased when n1 + n2 ≥ N, with its exact bias otherwise, by enumerating the hypergeometric law;
  - Schnabel by hand, and equal to Lincoln–Petersen on two occasions;
  - M0 and Mt equal exhaustive searches over N, and the likelihood maximised over p on a grid;
  - 300 simulated populations each: bias and the profile intervals' coverage; AIC between M0 and Mt.
- `check_peira_distance_sampling` (32):
  - the half-normal's closed forms without truncation (σ̂² = Σx²/n, Σr²/(2n); μ = σ√(π/2); cv 1/√(2n) and 1/√n);
  - the truncated closed forms against quadrature, and the hazard rate's quadrature against σΓ(1 − 1/b);
  - D̂ × covered area equals `weighting`'s Horvitz–Thompson total;
  - R2 by hand;
  - cosine adjustments;
  - 200 line surveys (half-normal) and 150 point surveys (hazard rate): bias, log-normal coverage and standard errors; AIC chooses the hazard rate on hazard-rate distances.

No reference package was run: R's `unmarked` and `Distance` are not installed. The checks rest on what can be computed independently of the code.

## 4. Dependencies

tyche `probability` (distributions, quantiles); calx `optimisation` (BFGS, Nelder–Mead, derivatives by differences), `special_functions` (log Γ, the incomplete beta function), `quadrature` (Gauss–Kronrod), `solvers`, `linalg` (Cholesky), `random` (the simulations' streams). The distance-sampling check also uses peira `weighting`.

## 5. References

- Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers, D. L. and Thomas, L. (2001), *Introduction to Distance Sampling*, Oxford University Press.
- Chao, A. (1987), Estimating the population size for capture–recapture data with unequal catchability, *Biometrics* 43, 783–791.
- Chapman, D. G. (1951), Some properties of the hypergeometric distribution with applications to zoological sample censuses, *Univ. California Publ. Statist.* 1, 131–160.
- Darroch, J. N. (1958), The multiple-recapture census I: estimation of a closed population, *Biometrika* 45, 343–359.
- Dennis, E. B., Morgan, B. J. T. and Ridout, M. S. (2015), Computational aspects of N-mixture models, *Biometrics* 71, 237–246.
- Fewster, R. M., Buckland, S. T., Burnham, K. P., Borchers, D. L., Jupp, P. E., Laake, J. L. and Thomas, L. (2009), Estimating the encounter rate variance in distance sampling, *Biometrics* 65, 225–236.
- Hayes, R. J. and Buckland, S. T. (1983), Radial-distance models for the line-transect method, *Biometrics* 39, 29–42.
- Lincoln, F. C. (1930), Calculating waterfowl abundance on the basis of banding returns, USDA Circular 118.
- Louis, T. A. (1982), Finding the observed information matrix when using the EM algorithm, *J. R. Statist. Soc. B* 44, 226–233.
- MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle, J. A. and Langtimm, C. A. (2002), Estimating site occupancy rates when detection probabilities are less than one, *Ecology* 83, 2248–2255.
- MacKenzie, D. I., Nichols, J. D., Royle, J. A., Pollock, K. H., Bailey, L. L. and Hines, J. E. (2018), *Occupancy Estimation and Modeling*, 2nd ed., Academic Press.
- Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978), Statistical inference from capture data on closed animal populations, *Wildlife Monographs* 62.
- Petersen, C. G. J. (1896), The yearly immigration of young plaice into the Limfjord from the German Sea, *Rep. Danish Biol. Stn* 6, 1–48.
- Royle, J. A. (2004), N-mixture models for estimating population size from spatially replicated counts, *Biometrics* 60, 108–115.
- Schnabel, Z. E. (1938), The estimation of the total fish population of a lake, *Amer. Math. Monthly* 45, 348–352.
- Seber, G. A. F. (1970), The effects of trap response on tag recapture estimates, *Biometrics* 26, 13–22.
- Seber, G. A. F. (1982), *The Estimation of Animal Abundance and Related Parameters*, 2nd ed., Griffin.
