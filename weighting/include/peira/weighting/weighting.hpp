// Weighting by the inverse of the probability of being seen (docs/theory.md, §2 and §4).
//
// One estimator under two names. A sample s drawn from a population U with inclusion probabilities π_i = P(i ∈ s) > 0
// reveals y_i for i ∈ s; Horvitz and Thompson (1952) estimate the total by Σ_{i∈s} y_i/π_i, unbiased for any design.
// A treatment assignment reveals Y_i(1) for the treated with probability e_i = P(Z_i = 1) and Y_i(0) for the controls with
// probability 1 − e_i; the same estimator of each arm's total gives the inverse-probability-weighted effect.
//
// - horvitz_thompson_total, hajek_mean: the total, and the ratio Σ y/π / Σ 1/π (Hájek 1971), which keeps the weights'
//   sum and trades a small bias for less variance.
// - ht_variance, syg_variance: the unbiased variance estimators of the total from the joint inclusion probabilities
//   π_ij > 0, in Horvitz and Thompson's form and, for designs of fixed size, Sen–Yates–Grundy's (Särndal, Swensson and
//   Wretman 1992, Model Assisted Survey Sampling, §2.8).
// - ipw: the effect by inverse-probability weights (normalised within each arm, or not); aipw: the augmented estimator of
//   Robins, Rotnitzky and Zhao (1994), with ψ_i = m̂₁ − m̂₀ + Z(Y − m̂₁)/ê − (1 − Z)(Y − m̂₀)/(1 − ê), consistent when the
//   propensity or the outcome model is right (doubly robust), its standard error from its influence function.
// - fit_nuisances: the propensity by logistic regression and each arm's outcome by least squares, through tyche.
//
// Limits: effects are identified only under unconfoundedness given the covariates and overlap (0 < e < 1); weights near
// zero or one make every weighted estimator unstable, and `trim` is offered, not applied by default.

#pragma once

#include <tyche/probability/distributions.hpp>
#include <tyche/regression/glm.hpp>
#include <tyche/regression/linear.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/rectangular.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace peira::weighting {

using calx::linalg::Matrix;
using calx::linalg::Rectangular;

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::weighting: ") + what);
}
inline double z_two_sided(double confidence) { return tyche::probability::Normal{0, 1}.upper_quantile((1 - confidence) / 2); }
}  // namespace detail

inline double horvitz_thompson_total(std::span<const double> y, std::span<const double> pi) {
  detail::require(y.size() == pi.size(), "one inclusion probability per unit");
  double t = 0;
  for (std::size_t i = 0; i < y.size(); i++) {
    detail::require(pi[i] > 0 && pi[i] <= 1, "inclusion probabilities in (0, 1]");
    t += y[i] / pi[i];
  }
  return t;
}

inline double hajek_mean(std::span<const double> y, std::span<const double> pi) {
  double t = 0, n = 0;
  for (std::size_t i = 0; i < y.size(); i++) { t += y[i] / pi[i]; n += 1 / pi[i]; }
  return t / n;
}

// Σ_i Σ_j (π_ij − π_i π_j)/π_ij · (y_i/π_i)(y_j/π_j), with π_ii = π_i: unbiased for Var of the total when every π_ij > 0.
inline double ht_variance(std::span<const double> y, std::span<const double> pi, const Matrix& pij) {
  const std::size_t n = y.size();
  detail::require(pij.n == n, "a joint inclusion probability for every pair in the sample");
  double v = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      const double p = i == j ? pi[i] : pij(i, j);
      v += (p - pi[i] * pi[j]) / p * (y[i] / pi[i]) * (y[j] / pi[j]);
    }
  return v;
}

// ½ Σ_{i≠j} (π_i π_j − π_ij)/π_ij · (y_i/π_i − y_j/π_j)²: unbiased for designs of fixed size.
inline double syg_variance(std::span<const double> y, std::span<const double> pi, const Matrix& pij) {
  const std::size_t n = y.size();
  detail::require(pij.n == n, "a joint inclusion probability for every pair in the sample");
  double v = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      if (i == j) continue;
      const double d = y[i] / pi[i] - y[j] / pi[j];
      v += 0.5 * (pi[i] * pi[j] - pij(i, j)) / pij(i, j) * d * d;
    }
  return v;
}

struct Effect {
  double estimate = 0, se = 0, ci_lower = 0, ci_upper = 0;
  std::vector<double> influence;  // ψ_i − τ̂, one per unit
};

namespace detail {
inline Effect from_influence(std::vector<double> psi, double confidence) {
  const double n = static_cast<double>(psi.size());
  Effect r;
  for (double v : psi) r.estimate += v / n;
  double ss = 0;
  for (double& v : psi) { v -= r.estimate; ss += v * v; }
  r.se = std::sqrt(ss / (n * (n - 1)));
  const double q = z_two_sided(confidence);
  r.ci_lower = r.estimate - q * r.se, r.ci_upper = r.estimate + q * r.se;
  r.influence = std::move(psi);
  return r;
}
}  // namespace detail

// The average effect by inverse-probability weights. normalised: each arm's weighted mean is Hájek's; else Horvitz and
// Thompson's (1/n) Σ [Z Y/e − (1 − Z) Y/(1 − e)]. The standard error treats e as known (its influence function).
inline Effect ipw(std::span<const double> y, std::span<const int> z, std::span<const double> e, bool normalised = true, double confidence = 0.95) {
  const std::size_t n = y.size();
  detail::require(z.size() == n && e.size() == n, "one assignment and one propensity per outcome");
  std::vector<double> psi(n);
  if (!normalised) {
    for (std::size_t i = 0; i < n; i++) {
      detail::require(e[i] > 0 && e[i] < 1, "propensities in (0, 1) (overlap)");
      psi[i] = z[i] ? y[i] / e[i] : -y[i] / (1 - e[i]);
    }
    return detail::from_influence(std::move(psi), confidence);
  }
  double s1 = 0, w1 = 0, s0 = 0, w0 = 0;
  for (std::size_t i = 0; i < n; i++) {
    detail::require(e[i] > 0 && e[i] < 1, "propensities in (0, 1) (overlap)");
    if (z[i]) { s1 += y[i] / e[i]; w1 += 1 / e[i]; } else { s0 += y[i] / (1 - e[i]); w0 += 1 / (1 - e[i]); }
  }
  const double mu1 = s1 / w1, mu0 = s0 / w0;
  // The linearised influence of the ratio estimators: n Z (Y − μ₁)/(e Σ Z/e) − n (1 − Z)(Y − μ₀)/((1 − e) Σ (1−Z)/(1−e)).
  const double nn = static_cast<double>(n);
  for (std::size_t i = 0; i < n; i++)
    psi[i] = (mu1 - mu0) + (z[i] ? nn * (y[i] - mu1) / (e[i] * w1) : -nn * (y[i] - mu0) / ((1 - e[i]) * w0));
  return detail::from_influence(std::move(psi), confidence);
}

// The augmented (doubly robust) estimator, given the propensities ê and the arms' predicted outcomes m̂₁, m̂₀.
inline Effect aipw(std::span<const double> y, std::span<const int> z, std::span<const double> e, std::span<const double> m1, std::span<const double> m0,
                   double confidence = 0.95) {
  const std::size_t n = y.size();
  detail::require(z.size() == n && e.size() == n && m1.size() == n && m0.size() == n, "one of each per unit");
  std::vector<double> psi(n);
  for (std::size_t i = 0; i < n; i++) {
    detail::require(e[i] > 0 && e[i] < 1, "propensities in (0, 1) (overlap)");
    psi[i] = m1[i] - m0[i] + (z[i] ? (y[i] - m1[i]) / e[i] : -(y[i] - m0[i]) / (1 - e[i]));
  }
  return detail::from_influence(std::move(psi), confidence);
}

// Propensities clipped to [a, 1 − a].
inline std::vector<double> trim(std::span<const double> e, double a = 0.01) {
  std::vector<double> out(e.begin(), e.end());
  for (double& v : out) v = std::clamp(v, a, 1 - a);
  return out;
}

struct Nuisances {
  std::vector<double> e, m1, m0;
};

// The propensity e(x) by a logistic regression of Z on P, and each arm's outcome by least squares of Y on Q within the
// arm, predicted for every unit. P and Q are design matrices with their intercepts.
inline Nuisances fit_nuisances(const Rectangular& P, const Rectangular& Q, std::span<const double> y, std::span<const int> z) {
  const std::size_t n = y.size();
  detail::require(P.rows == n && Q.rows == n && z.size() == n, "one row of each design matrix per unit");
  std::vector<double> zd(n);
  for (std::size_t i = 0; i < n; i++) zd[i] = z[i];
  tyche::regression::GlmOptions o;
  o.family = tyche::regression::Family::binomial;
  const auto g = tyche::regression::glm(P, zd, o);
  Nuisances r;
  r.e = g.fitted;
  for (int arm = 0; arm <= 1; arm++) {
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; i++) k += z[i] == arm;
    Rectangular Qa(k, Q.cols);
    std::vector<double> ya(k);
    std::size_t row = 0;
    for (std::size_t i = 0; i < n; i++) {
      if (z[i] != arm) continue;
      for (std::size_t j = 0; j < Q.cols; j++) Qa(row, j) = Q(i, j);
      ya[row++] = y[i];
    }
    const auto fit = tyche::regression::ols(Qa, ya);
    std::vector<double>& m = arm ? r.m1 : r.m0;
    m.assign(n, 0.0);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < Q.cols; j++) m[i] += Q(i, j) * fit.beta[j];
  }
  return r;
}

}  // namespace peira::weighting
