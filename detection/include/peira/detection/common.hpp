// What the designs of imperfect detection share (docs/theory.md, §10; peira's decision 0008).
//
// A survey of sites on repeated visits records, at site i on visit j, a count y_ij (a detection, 0 or 1, when what is recorded
// is presence) or nothing, when the site was not visited then. SiteVisits holds that table, with −1 for "not visited". The
// models fitted to it are logit-linear or log-linear in covariates given as design matrices: one row per site for what the site
// holds (its occupancy ψ_i or its abundance λ_i), one row per site and visit, row i·visits + j, for the probability p_ij of
// detecting what is there. Design matrices carry their own intercepts; `intercept` makes the constant one.
//
// - predict: a fitted probability or mean at a row of covariates, x'θ through the inverse link, with the delta method's
//   standard error (the gradient of the inverse link times x, against the coefficients' covariance) and the interval taken on
//   the link's scale and carried through it, so that it stays within (0, 1) or (0, ∞).
// - detail::information_inverse: the inverse of the observed information −∇²ℓ, by Cholesky, or nothing when it is not positive
//   definite (an estimate on the boundary, or a model not identified by the data).

#pragma once

#include <tyche/probability/distributions.hpp>

#include <calx/linalg/cholesky.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/linalg/rectangular.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace peira::detection {

using calx::linalg::Matrix;
using calx::linalg::Rectangular;

inline constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::detection: ") + what);
}
inline double z_two_sided(double confidence) { return tyche::probability::Normal{0, 1}.upper_quantile((1 - confidence) / 2); }

// log(1 + eˣ), without overflow for large x or loss for very negative x.
inline double log1pexp(double x) { return x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x)); }
// 1/(1 + e⁻ˣ), without overflow.
inline double logistic(double x) {
  if (x >= 0) return 1 / (1 + std::exp(-x));
  const double e = std::exp(x);
  return e / (1 + e);
}

// (−H)⁻¹ when −H is positive definite; false otherwise.
inline bool information_inverse(const Matrix& hessian, Matrix& out) {
  const std::size_t n = hessian.n;
  Matrix info(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) info(i, j) = -hessian(i, j);
  const auto f = calx::linalg::cholesky(info);
  if (!f.positive_definite) return false;
  out = Matrix(n);
  std::vector<double> e(n, 0.0);
  for (std::size_t j = 0; j < n; j++) {
    std::fill(e.begin(), e.end(), 0.0);
    e[j] = 1;
    const auto col = calx::linalg::cholesky_solve(f, e);
    for (std::size_t i = 0; i < n; i++) out(i, j) = col[i];
  }
  return true;
}

inline std::vector<double> standard_errors(const Matrix& covariance, bool available, std::size_t n) {
  std::vector<double> se(n, not_a_number);
  if (available)
    for (std::size_t i = 0; i < n; i++) se[i] = std::sqrt(covariance(i, i));
  return se;
}
}  // namespace detail

// What was seen at sites on visits; y(i, j) = −1 where site i was not visited on visit j.
struct SiteVisits {
  std::size_t sites = 0, visits = 0;
  std::vector<int> y;  // row-major: y[i * visits + j]

  SiteVisits() = default;
  SiteVisits(std::size_t s, std::size_t v) : sites(s), visits(v), y(s * v, -1) {}
  int& operator()(std::size_t i, std::size_t j) { return y[i * visits + j]; }
  int operator()(std::size_t i, std::size_t j) const { return y[i * visits + j]; }
};

// A design matrix of one column of ones: the constant model.
inline Rectangular intercept(std::size_t rows) {
  Rectangular m(rows, 1);
  for (std::size_t i = 0; i < rows; i++) m(i, 0) = 1;
  return m;
}

// How a model is fitted by maximum likelihood.
struct FitOptions {
  double gradient_tolerance = 1e-9;  // on the mean log-likelihood per site (or per detection, for distances), in the max norm
  int max_iterations = 0;            // BFGS's; 0: its default
  std::vector<double> start;         // the coefficients to start from; empty: each model's default
};

enum class Link { logit, log };

struct Prediction {
  double estimate = 0, se = 0, lower = 0, upper = 0;  // on the scale of the probability or the mean
  double linear = 0, linear_se = 0;                    // x'θ and its standard error
};

// The inverse link at x'θ, θ the coefficients at positions offset … offset + x.size() − 1 of a fit whose covariance is given.
inline Prediction predict(std::span<const double> theta, const Matrix& covariance, std::size_t offset, std::span<const double> x, Link link,
                          double confidence = 0.95) {
  detail::require(offset + x.size() <= theta.size() && covariance.n == theta.size(), "a row of covariates within the fit's coefficients");
  Prediction r;
  for (std::size_t a = 0; a < x.size(); a++) r.linear += x[a] * theta[offset + a];
  double v = 0;
  for (std::size_t a = 0; a < x.size(); a++)
    for (std::size_t b = 0; b < x.size(); b++) v += x[a] * x[b] * covariance(offset + a, offset + b);
  r.linear_se = std::sqrt(std::max(0.0, v));
  const double q = detail::z_two_sided(confidence);
  auto inverse = [&](double eta) { return link == Link::logit ? detail::logistic(eta) : std::exp(eta); };
  r.estimate = inverse(r.linear);
  const double slope = link == Link::logit ? r.estimate * (1 - r.estimate) : r.estimate;
  r.se = slope * r.linear_se;
  r.lower = inverse(r.linear - q * r.linear_se);
  r.upper = inverse(r.linear + q * r.linear_se);
  return r;
}

}  // namespace peira::detection
