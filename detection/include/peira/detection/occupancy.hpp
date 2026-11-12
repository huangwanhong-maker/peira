// Site occupancy with imperfect detection (docs/theory.md, §10.2): MacKenzie, Nichols, Lachman, Droege, Royle and Langtimm
// (2002), "Estimating site occupancy rates when detection probabilities are less than one", Ecology 83, 2248–2255; MacKenzie,
// Nichols, Royle, Pollock, Bailey and Hines (2018), Occupancy Estimation and Modeling, 2nd ed., Academic Press (the single-season
// model).
//
// FULL DATA AND MECHANISM. Each site i is occupied (z_i = 1) or not (z_i = 0), with P(z_i = 1) = ψ_i = logistic(x_i'α). The
// occupancy is the full data; the mechanism is detection: on visit j an occupied site is detected with probability
// p_ij = logistic(w_ij'β), an unoccupied one never. What is seen is the detection history y_i. The likelihood is
//   L_i = ψ_i Π_j p_ij^{y_ij} (1 − p_ij)^{1 − y_ij} + (1 − ψ_i) [y_i = 0],
// the case N ∈ {0, 1} of latent.hpp. A site never detected is either empty or missed, and the mixture separates the two only
// through replication: with one visit the likelihood depends on (ψ, p) only through ψ p, and ψ is not identified. With two or
// more visits and ψ, p constant, the histories of the detected sites follow a binomial truncated at zero, which identifies p,
// and the share of sites detected then gives ψ.
//
// IDENTIFICATION STATED. (i) Closure: z_i does not change between visits (untestable within the design: a site that empties and
// fills reads as one detected imperfectly). (ii) No false detections. (iii) Detections independent given z_i and the covariates
// (untestable; unmodelled heterogeneity in p biases ψ̂ low). (iv) The models for ψ and p are right.
//
// - occupancy: the maximum likelihood estimates, with standard errors from the observed information (Louis's identity, exact);
//   AIC = −2ℓ + 2k; the fitted ψ_i and p_ij; the posterior P(z_i = 1 | y_i), 1 where the site was detected and otherwise
//   ψ_i Π_j (1 − p_ij) / (ψ_i Π_j (1 − p_ij) + 1 − ψ_i), by Bayes's rule; the naive occupancy, the share of
//   sites surveyed with a detection, which estimates ψ p* and not ψ; and the finite-sample occupancy, the expected share of the
//   sites surveyed that are occupied given what was seen, Σ_i P(z_i = 1 | y_i) / s (an estimate of the realised occupancy of
//   these sites, not of ψ of their population).
// - predict_psi, predict_p: a probability at a row of covariates, with the delta method's standard error.
//
// Limits: estimates near ψ = 1 or p = 0 lie on the boundary, where the information is singular and the standard errors are not
// available (information_positive = false); the asymptotics need many sites, not many visits.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include "common.hpp"
#include "latent.hpp"

namespace peira::detection {

struct Occupancy {
  std::vector<double> psi_coefficients, psi_se;  // α: ψ_i = logistic(x_i'α)
  std::vector<double> p_coefficients, p_se;      // β: p_ij = logistic(w_ij'β)
  std::vector<double> theta;                     // (α, β)
  Matrix covariance;                             // of θ, the inverse observed information
  double log_likelihood = 0, aic = 0;
  std::size_t parameters = 0;
  bool converged = false, information_positive = false;
  int iterations = 0;
  std::vector<double> psi;       // ψ_i, per site
  std::vector<double> p;         // p_ij, per site and visit (row-major); NaN where the visit was not made
  std::vector<double> occupied;  // P(z_i = 1 | y_i), per site; ψ_i where the site was never visited
  std::size_t sites_surveyed = 0, sites_detected = 0;
  double naive = 0;          // sites detected / sites surveyed
  double finite_sample = 0;  // Σ P(z_i = 1 | y_i) over the sites surveyed / sites surveyed
  double mean_psi = 0;       // the mean of ψ_i over the sites surveyed
};

// ℓ at θ = (α, β), for profiles and checks.
inline double occupancy_log_likelihood(const SiteVisits& y, const Rectangular& X, const Rectangular& W, std::span<const double> theta) {
  detail::latent_validate(y, X, W, Mixing::bernoulli, 1);
  return detail::latent_evaluate(y, X, W, Mixing::bernoulli, 1, theta, 0).loglik;
}

// X: one row per site (ψ's covariates, with an intercept); W: one row per site and visit (p's covariates), row i·visits + j.
// The default start: ψ at the naive occupancy (on X's intercept, if it has one), p at 1/2.
inline Occupancy occupancy(const SiteVisits& y, const Rectangular& X, const Rectangular& W, const FitOptions& o = {}) {
  detail::latent_validate(y, X, W, Mixing::bernoulli, 1);
  const std::size_t S = y.sites, T = y.visits, kx = X.cols, kw = W.cols;
  Occupancy r;
  for (std::size_t i = 0; i < S; i++) {
    bool visited = false, detected = false;
    for (std::size_t j = 0; j < T; j++) visited |= y(i, j) >= 0, detected |= y(i, j) == 1;
    r.sites_surveyed += visited, r.sites_detected += detected;
  }
  detail::require(r.sites_surveyed > 0, "at least one site visited");
  r.naive = static_cast<double>(r.sites_detected) / static_cast<double>(r.sites_surveyed);
  std::vector<double> start = o.start;
  if (start.empty()) {
    start.assign(kx + kw, 0.0);
    const double v = std::clamp(r.naive, 0.05, 0.95);
    detail::start_on_constant_columns(X, std::log(v / (1 - v)), start, 0);
  }
  const auto fit = detail::latent_fit(y, X, W, Mixing::bernoulli, 1, std::move(start), o);
  r.theta = fit.theta;
  r.covariance = fit.covariance;
  r.psi_coefficients.assign(fit.theta.begin(), fit.theta.begin() + static_cast<std::ptrdiff_t>(kx));
  r.p_coefficients.assign(fit.theta.begin() + static_cast<std::ptrdiff_t>(kx), fit.theta.end());
  r.psi_se.assign(fit.se.begin(), fit.se.begin() + static_cast<std::ptrdiff_t>(kx));
  r.p_se.assign(fit.se.begin() + static_cast<std::ptrdiff_t>(kx), fit.se.end());
  r.log_likelihood = fit.loglik;
  r.parameters = kx + kw;
  r.aic = -2 * r.log_likelihood + 2 * static_cast<double>(r.parameters);
  r.converged = fit.converged;
  r.information_positive = fit.information_positive;
  r.iterations = fit.iterations;
  std::vector<std::vector<double>> post;
  detail::latent_evaluate(y, X, W, Mixing::bernoulli, 1, fit.theta, 0, &post);
  r.psi.resize(S), r.occupied.resize(S), r.p.assign(S * T, not_a_number);
  double sum_post = 0, sum_psi = 0;
  for (std::size_t i = 0; i < S; i++) {
    double eta = 0;
    for (std::size_t a = 0; a < kx; a++) eta += X(i, a) * fit.theta[a];
    r.psi[i] = detail::logistic(eta);
    r.occupied[i] = post[i][1];
    bool visited = false;
    for (std::size_t j = 0; j < T; j++) {
      if (y(i, j) < 0) continue;
      visited = true;
      double ed = 0;
      for (std::size_t b = 0; b < kw; b++) ed += W(i * T + j, b) * fit.theta[kx + b];
      r.p[i * T + j] = detail::logistic(ed);
    }
    if (visited) sum_post += r.occupied[i], sum_psi += r.psi[i];
  }
  r.finite_sample = sum_post / static_cast<double>(r.sites_surveyed);
  r.mean_psi = sum_psi / static_cast<double>(r.sites_surveyed);
  return r;
}

inline Prediction predict_psi(const Occupancy& fit, std::span<const double> x, double confidence = 0.95) {
  detail::require(x.size() == fit.psi_coefficients.size(), "one covariate per coefficient of ψ");
  return predict(fit.theta, fit.covariance, 0, x, Link::logit, confidence);
}
inline Prediction predict_p(const Occupancy& fit, std::span<const double> w, double confidence = 0.95) {
  detail::require(w.size() == fit.p_coefficients.size(), "one covariate per coefficient of p");
  return predict(fit.theta, fit.covariance, fit.psi_coefficients.size(), w, Link::logit, confidence);
}

}  // namespace peira::detection
