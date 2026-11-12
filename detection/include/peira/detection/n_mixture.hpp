// Abundance from repeated counts: the N-mixture model (docs/theory.md, §10.3). Royle (2004), "N-mixture models for estimating
// population size from spatially replicated counts", Biometrics 60, 108–115.
//
// FULL DATA AND MECHANISM. Site i holds N_i individuals, N_i ~ Poisson(λ_i) or negative binomial with mean λ_i and size r,
// λ_i = exp(x_i'α); that is the full data. On visit j each is counted independently with probability p_ij = logistic(w_ij'β), so
// the count is y_ij ~ Binomial(N_i, p_ij); what is seen is the counts. The likelihood sums over N_i ≥ max_j y_ij (latent.hpp),
// truncated at a bound K that the caller states. With one visit the Poisson model is not identified: Σ_N Poisson(N; λ)
// Binomial(y; N, p) = Poisson(y; λp), a function of λp alone (thinning). Repeated counts at a site, with N_i fixed between
// them, identify p from their agreement.
//
// IDENTIFICATION STATED. (i) Closure: N_i does not change between visits (untestable within the design; movement in and out
// reads as a larger λ and a smaller p). (ii) No double counts and no false counts. (iii) Individuals detected independently
// with a common p_ij (untestable: heterogeneity in p biases λ̂). (iv) The mixing law is right: the Poisson and the negative
// binomial can give very different λ̂ from the same counts, and with few visits the likelihood can be maximised at λ = ∞
// with positive probability (Dennis, Morgan and Ridout 2015, Biometrics 71, 237–246), so the estimates' stability under the
// choice of law and of K is part of every analysis.
//
// THE BOUND K. The sum over N stops at K. The fit reports two measures of what is lost: the largest prior probability beyond
// K at any site under the fitted law, max_i P(N_i > K) (Poisson by tyche; negative binomial as I_{λ/(r+λ)}(K + 1, r) by calx's
// incomplete beta function), and the change in the log-likelihood at the estimate when K is doubled. Both should be negligible;
// refitting with a larger K is the check.
//
// - n_mixture: the maximum likelihood estimates, with standard errors from the observed information (exact, by Louis's
//   identity); AIC = −2ℓ + 2k; the fitted λ_i and p_ij; the posterior law of each N_i given its counts, N = 0 … K, and its mean,
//   the best predictor of the site's abundance; the totals Σ λ_i and Σ E[N_i | y_i] over the sites surveyed.
// - predict_lambda, predict_p: a mean or a probability at a row of covariates, with the delta method's standard error.

#pragma once

#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include "common.hpp"
#include "latent.hpp"

namespace peira::detection {

struct NMixture {
  Mixing mixing = Mixing::poisson;
  int K = 0;
  std::vector<double> lambda_coefficients, lambda_se;  // α: λ_i = exp(x_i'α)
  std::vector<double> p_coefficients, p_se;            // β: p_ij = logistic(w_ij'β)
  double log_size = not_a_number, log_size_se = not_a_number;  // ρ = log r (negative binomial only)
  double size = not_a_number;                                   // r
  std::vector<double> theta;                           // (α, β[, ρ])
  Matrix covariance;
  double log_likelihood = 0, aic = 0;
  std::size_t parameters = 0;
  bool converged = false, information_positive = false;
  int iterations = 0;
  std::vector<double> lambda;                    // λ_i, per site
  std::vector<double> p;                         // p_ij, row-major; NaN where the visit was not made
  std::vector<std::vector<double>> posterior;    // per site, P(N_i = N | y_i) for N = 0 … K (the prior, truncated, where never visited)
  std::vector<double> abundance;                 // E[N_i | y_i]
  std::size_t sites_surveyed = 0;
  double total_expected = 0, total_predicted = 0;  // Σ λ_i and Σ E[N_i | y_i] over the sites surveyed
  double tail_probability = 0;                     // max over the sites surveyed of P(N_i > K) under the fitted law
  double truncation_effect = 0;                    // ℓ(θ̂; 2K) − ℓ(θ̂; K)
};

// ℓ at θ = (α, β[, ρ]) with the sum truncated at K.
inline double n_mixture_log_likelihood(const SiteVisits& y, const Rectangular& X, const Rectangular& W, int K, Mixing mixing, std::span<const double> theta) {
  detail::require(mixing != Mixing::bernoulli, "a law of abundance (Poisson or negative binomial)");
  detail::latent_validate(y, X, W, mixing, K);
  return detail::latent_evaluate(y, X, W, mixing, K, theta, 0).loglik;
}

// X: one row per site (λ's covariates, with an intercept); W: one row per site and visit (p's), row i·visits + j; K ≥ every count.
// The default start: λ at twice the mean of the sites' largest counts (on X's intercept, if it has one), p at 1/2, r at 5.
inline NMixture n_mixture(const SiteVisits& y, const Rectangular& X, const Rectangular& W, int K, Mixing mixing = Mixing::poisson, const FitOptions& o = {}) {
  detail::require(mixing != Mixing::bernoulli, "a law of abundance (Poisson or negative binomial); presence is occupancy()");
  detail::latent_validate(y, X, W, mixing, K);
  const std::size_t S = y.sites, T = y.visits, kx = X.cols, kw = W.cols, k = detail::latent_parameters(X, W, mixing);
  const bool nb = mixing == Mixing::negative_binomial;
  NMixture r;
  r.mixing = mixing, r.K = K;
  std::vector<bool> visited(S, false);
  double sum_max = 0;
  for (std::size_t i = 0; i < S; i++) {
    int m = -1;
    for (std::size_t j = 0; j < T; j++) m = std::max(m, y(i, j));
    if (m >= 0) visited[i] = true, r.sites_surveyed++, sum_max += m;
  }
  detail::require(r.sites_surveyed > 0, "at least one site visited");
  std::vector<double> start = o.start;
  if (start.empty()) {
    start.assign(k, 0.0);
    detail::start_on_constant_columns(X, std::log(std::max(0.5, 2 * sum_max / static_cast<double>(r.sites_surveyed))), start, 0);
    if (nb) start[k - 1] = std::log(5.0);
  }
  const auto fit = detail::latent_fit(y, X, W, mixing, K, std::move(start), o);
  r.theta = fit.theta;
  r.covariance = fit.covariance;
  r.lambda_coefficients.assign(fit.theta.begin(), fit.theta.begin() + static_cast<std::ptrdiff_t>(kx));
  r.p_coefficients.assign(fit.theta.begin() + static_cast<std::ptrdiff_t>(kx), fit.theta.begin() + static_cast<std::ptrdiff_t>(kx + kw));
  r.lambda_se.assign(fit.se.begin(), fit.se.begin() + static_cast<std::ptrdiff_t>(kx));
  r.p_se.assign(fit.se.begin() + static_cast<std::ptrdiff_t>(kx), fit.se.begin() + static_cast<std::ptrdiff_t>(kx + kw));
  if (nb) r.log_size = fit.theta[k - 1], r.log_size_se = fit.se[k - 1], r.size = std::exp(r.log_size);
  r.log_likelihood = fit.loglik;
  r.parameters = k;
  r.aic = -2 * r.log_likelihood + 2 * static_cast<double>(k);
  r.converged = fit.converged;
  r.information_positive = fit.information_positive;
  r.iterations = fit.iterations;
  detail::latent_evaluate(y, X, W, mixing, K, fit.theta, 0, &r.posterior);
  r.truncation_effect = detail::latent_evaluate(y, X, W, mixing, 2 * K, fit.theta, 0).loglik - fit.loglik;
  r.lambda.resize(S), r.abundance.resize(S), r.p.assign(S * T, not_a_number);
  for (std::size_t i = 0; i < S; i++) {
    double eta = 0;
    for (std::size_t a = 0; a < kx; a++) eta += X(i, a) * fit.theta[a];
    r.lambda[i] = std::exp(eta);
    double m = 0;
    for (std::size_t n = 0; n < r.posterior[i].size(); n++) m += static_cast<double>(n) * r.posterior[i][n];
    r.abundance[i] = m;
    for (std::size_t j = 0; j < T; j++) {
      if (y(i, j) < 0) continue;
      double ed = 0;
      for (std::size_t b = 0; b < kw; b++) ed += W(i * T + j, b) * fit.theta[kx + b];
      r.p[i * T + j] = detail::logistic(ed);
    }
    if (!visited[i]) continue;
    r.total_expected += r.lambda[i], r.total_predicted += m;
    const double tail = !std::isfinite(r.lambda[i]) || !(r.lambda[i] > 0) ? (r.lambda[i] > 0 ? 1.0 : 0.0)
                        : nb ? calx::special_functions::beta_i(static_cast<double>(K) + 1, r.size, r.lambda[i] / (r.size + r.lambda[i]))
                             : tyche::probability::Poisson{r.lambda[i]}.survival(K);
    r.tail_probability = std::max(r.tail_probability, tail);
  }
  return r;
}

inline Prediction predict_lambda(const NMixture& fit, std::span<const double> x, double confidence = 0.95) {
  detail::require(x.size() == fit.lambda_coefficients.size(), "one covariate per coefficient of λ");
  return predict(fit.theta, fit.covariance, 0, x, Link::log, confidence);
}
inline Prediction predict_p(const NMixture& fit, std::span<const double> w, double confidence = 0.95) {
  detail::require(w.size() == fit.p_coefficients.size(), "one covariate per coefficient of p");
  return predict(fit.theta, fit.covariance, fit.lambda_coefficients.size(), w, Link::logit, confidence);
}

}  // namespace peira::detection
