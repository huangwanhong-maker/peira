// Distance sampling (docs/theory.md, §10.5): Buckland, Anderson, Burnham, Laake, Borchers and Thomas (2001), Introduction to
// Distance Sampling, Oxford University Press, ch. 3; the variance of the encounter rate after Fewster, Buckland, Burnham,
// Borchers, Jupp, Laake and Thomas (2009), "Estimating the encounter rate variance in distance sampling", Biometrics 65, 225–236.
//
// FULL DATA AND MECHANISM. Individuals lie in the plane with density D; that is the full data. An observer walks lines of total
// length L (or stands at k points) and detects an individual at perpendicular (or radial) distance x with probability g(x),
// the detection function, with g(0) = 1. What is seen is the n individuals detected within the truncation distance w and
// their distances. Because the lines or points are placed at random, the individuals within w lie uniformly in distance from
// the line (density proportional to r from a point), so the distances seen have the density
//   lines:  f(x) = g(x)/μ,          μ = ∫₀^w g(x) dx, the effective strip half-width;
//   points: f(r) = 2π r g(r)/ν,     ν = 2π ∫₀^w r g(r) dr, the effective area of a point,
// and the probability that an individual in the covered region is detected is P_a = μ/w (lines) or ν/(πw²) (points). Then
//   D̂ = n f̂(0)/(2L) = n/(2 L μ̂)   (lines),      D̂ = n/(k ν̂)   (points),
// and both are Horvitz–Thompson estimates, n/P̂_a individuals in the covered area 2wL or kπw² (§10.1 of the theory note).
//
// IDENTIFICATION STATED. (i) g(0) = 1: everything on the line or at the point is seen (untestable from the distances: they give
// the shape of g, not its height). (ii) Lines or points placed independently of the individuals, so that their distances are
// uniform (or triangular) in the covered region; with a design that ensures it this holds by construction. (iii) Distances
// measured without error, and individuals detected before they move. (iv) The detection function's form.
//
// DETECTION FUNCTIONS (Buckland et al. 2001, ch. 3), a key with optional cosine adjustments:
// - uniform: g = 1 (with adjustments, the Fourier series estimator);
// - half_normal: g(x) = exp(−x²/(2σ²));
// - hazard_rate: g(x) = 1 − exp(−(x/σ)^{−b}) (Hayes and Buckland 1983), with a shoulder at 0;
// - cosine adjustments of orders j: g(x) = key(x) (1 + Σ_j a_j cos(jπx/w)) / (1 + Σ_j a_j), so that g(0) = 1; a fit that
//   makes g negative anywhere on [0, w] (on a grid of 129 points and at every distance) is refused. Monotonicity is not
//   imposed.
// The parameters are θ = (log σ [, log b], a_j …). μ and ν are closed forms for the half-normal without adjustments
// (σ √(π/2) erf(w/(σ√2)) and 2πσ²(1 − e^{−w²/(2σ²)})), and otherwise calx's adaptive Gauss–Kronrod quadrature, to a relative
// 10⁻¹². Without truncation (w = ∞) only the half-normal key without adjustments is allowed, whose MLE is then
// σ̂² = Σx²/n (lines) and Σr²/(2n) (points), exactly.
//
// ESTIMATION. Maximum likelihood by calx's Nelder–Mead from the start (three starts in b for the hazard rate), then its BFGS,
// then Newton's steps on calx's gradient and Hessian by central differences while the gradient falls, so that the estimate is
// found to about ε^{2/3} rather than the √ε at which searches on values stop. The covariance is the inverse of the Hessian by
// differences (relative error about 10⁻⁷ on smooth likelihoods).
// AIC = −2ℓ + 2q chooses among keys and adjustments. The cv of μ̂ (or ν̂) is the delta method's, ∇μ'Σ∇μ, with ∇μ by differences.
//
// DENSITY AND ITS VARIANCE. With n_k detections within w on transect k of effort e_k (a length for lines; for points, the number
// of visits, each visit counting once), E = Σ e_k and K transects:
//   encounter rate n/E, its variance (Fewster et al. 2009, the estimator R2)  K/(E²(K − 1)) Σ_k e_k² (n_k/e_k − n/E)²,
// applied in the same form to points with visits as effort; cv(D̂)² = cv(n/E)² + cv(μ̂)² (the delta method, encounter rate and
// detection function independent; Buckland et al. 2001, ch. 3); the log-normal interval D̂/C … D̂ C, C = exp(t √log(1 + cv²)),
// with t on Satterthwaite's degrees of freedom cv⁴/(cv_er⁴/(K − 1) + cv_det⁴/(n − q)) (Buckland et al. 2001, ch. 3).
// Distances and efforts are in one unit of length; D̂ is per that unit squared (individuals per m² when they are metres).
//
// Not written: grouped distances, clusters and their sizes, g(0) < 1 (double observers), covariates in the scale, monotonicity
// constraints on adjustments, stratification.

#pragma once

#include <calx/optimisation/nelder_mead.hpp>
#include <calx/optimisation/quasi_newton.hpp>
#include <calx/quadrature/gauss_kronrod.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include "common.hpp"

namespace peira::detection {

enum class Transect { line, point };
enum class Key { uniform, half_normal, hazard_rate };

struct DetectionModel {
  Key key = Key::half_normal;
  std::vector<int> cosine;  // the orders j ≥ 1 of the cosine adjustment terms; empty: the key alone
};

struct DistanceSurvey {
  Transect transect = Transect::line;
  std::vector<double> distance;   // one per detection: perpendicular (lines) or radial (points), ≥ 0
  std::vector<std::size_t> on;    // the transect each detection was made from, 0 … effort.size() − 1
  std::vector<double> effort;     // per transect: its length (lines), or the number of visits to it (points)
  double truncation = std::numeric_limits<double>::infinity();  // w; detections beyond it are left out
};

struct DetectionFit {
  DetectionModel model;
  Transect transect = Transect::line;
  double truncation = 0;
  std::vector<double> theta, se;  // (log σ [, log b], a_j …)
  Matrix covariance;
  double log_likelihood = 0, aic = 0;
  std::size_t n = 0, parameters = 0;
  bool converged = false, information_positive = false;
  double sigma = not_a_number, shape = not_a_number;  // σ and b, where the key has them
  double effective = 0;     // lines: μ, the effective strip half-width; points: ν, the effective area
  double effective_se = 0, effective_cv = 0;
  double f0 = 0;            // lines: f(0) = 1/μ; points: h(0) = 2π/ν
  double probability = not_a_number;  // P_a = μ/w or ν/(πw²); NaN without truncation
};

namespace detail {

inline std::size_t key_parameters(Key k) { return k == Key::uniform ? 0 : k == Key::half_normal ? 1 : 2; }

inline double key_g(Key key, std::span<const double> theta, double x) {
  switch (key) {
    case Key::uniform: return 1;
    case Key::half_normal: {
      const double s = std::exp(theta[0]);
      return std::exp(-x * x / (2 * s * s));
    }
    default: {
      const double s = std::exp(theta[0]), b = std::exp(theta[1]);
      return -std::expm1(-std::pow(x / s, -b));
    }
  }
}

inline double detection_g(const DetectionModel& m, std::span<const double> theta, double x, double w) {
  double g = key_g(m.key, theta, x);
  if (!m.cosine.empty()) {
    const std::size_t q = key_parameters(m.key);
    double s = 1, s0 = 1;
    for (std::size_t a = 0; a < m.cosine.size(); a++) {
      s += theta[q + a] * std::cos(m.cosine[a] * std::numbers::pi * x / w);
      s0 += theta[q + a];
    }
    g *= s / s0;
  }
  return g;
}

// μ (lines) or ν (points).
inline double effective(const DetectionModel& m, Transect t, std::span<const double> theta, double w) {
  if (m.key == Key::half_normal && m.cosine.empty()) {
    const double s = std::exp(theta[0]);
    if (t == Transect::line) return std::isfinite(w) ? s * std::sqrt(std::numbers::pi / 2) * std::erf(w / (s * std::numbers::sqrt2)) : s * std::sqrt(std::numbers::pi / 2);
    return 2 * std::numbers::pi * s * s * (std::isfinite(w) ? -std::expm1(-w * w / (2 * s * s)) : 1.0);
  }
  auto integrand = [&](double x) { return (t == Transect::line ? 1.0 : x) * detection_g(m, theta, x, w); };
  const auto e = calx::quadrature::integrate(integrand, 0.0, w, 1e-300, 1e-12, 2000);
  return (t == Transect::line ? 1.0 : 2 * std::numbers::pi) * e.value;
}

inline double distance_log_likelihood(const DetectionModel& m, Transect t, std::span<const double> x, double w, std::span<const double> theta) {
  const std::size_t q = key_parameters(m.key);
  for (std::size_t a = 0; a < q; a++)
    if (!(std::fabs(theta[a]) < 50)) return -std::numeric_limits<double>::infinity();  // σ and b within e^{±50}
  if (!m.cosine.empty()) {
    double s0 = 1;
    for (std::size_t a = 0; a < m.cosine.size(); a++) s0 += theta[q + a];
    if (!(s0 > 0)) return -std::numeric_limits<double>::infinity();
    for (int k = 0; k <= 128; k++)
      if (!(detection_g(m, theta, w * k / 128.0, w) >= 0)) return -std::numeric_limits<double>::infinity();
  }
  const double mu = effective(m, t, theta, w);
  if (!(mu > 0) || !std::isfinite(mu)) return -std::numeric_limits<double>::infinity();
  double l = -static_cast<double>(x.size()) * std::log(mu);
  for (double v : x) {
    const double g = detection_g(m, theta, v, w);
    if (!(g > 0)) return -std::numeric_limits<double>::infinity();
    l += std::log(g);
    if (t == Transect::point && v > 0) l += std::log(2 * std::numbers::pi * v);
  }
  return l;
}

inline std::vector<double> within(const DistanceSurvey& s) {
  std::vector<double> x;
  for (double v : s.distance) {
    require(v >= 0, "distances ≥ 0");
    if (v <= s.truncation) x.push_back(v);
  }
  return x;
}

}  // namespace detail

// g(x) of a fitted detection function.
inline double detection_probability(const DetectionFit& fit, double x) { return detail::detection_g(fit.model, fit.theta, x, fit.truncation); }

inline DetectionFit fit_detection(const DistanceSurvey& survey, const DetectionModel& model, const FitOptions& o = {}) {
  const double w = survey.truncation;
  detail::require(w > 0, "a truncation distance w > 0 (or ∞)");
  detail::require(std::isfinite(w) || (model.key == Key::half_normal && model.cosine.empty()),
                  "without truncation, only the half-normal key without adjustments");
  for (int j : model.cosine) detail::require(j >= 1, "cosine adjustments of orders ≥ 1");
  const auto x = detail::within(survey);
  const std::size_t n = x.size(), q = detail::key_parameters(model.key), k = q + model.cosine.size();
  detail::require(n > k, "more detections within w than parameters");
  DetectionFit r;
  r.model = model, r.transect = survey.transect, r.truncation = w, r.n = n, r.parameters = k;
  auto loglik = [&](std::span<const double> t) { return detail::distance_log_likelihood(model, survey.transect, x, w, t); };
  if (k == 0) {  // the uniform key alone: nothing to estimate
    r.log_likelihood = loglik({});
    r.converged = r.information_positive = true;
  } else {
    double m2 = 0;
    for (double v : x) m2 += v * v;
    m2 /= static_cast<double>(n);
    const double s0 = std::sqrt(survey.transect == Transect::line ? m2 : m2 / 2);
    std::vector<std::vector<double>> starts;
    if (!o.start.empty()) starts.push_back(o.start);
    else if (model.key == Key::hazard_rate)
      for (double b : {1.5, 3.0, 6.0}) {
        std::vector<double> t(k, 0.0);
        t[0] = std::log(0.8 * s0), t[1] = std::log(b);
        starts.push_back(t);
      }
    else {
      std::vector<double> t(k, 0.0);
      if (q == 1) t[0] = std::log(s0);
      starts.push_back(t);
    }
    const double scale = 1.0 / static_cast<double>(n);
    auto f = [&](std::span<const double> t) {
      const double l = loglik(t);
      return std::isfinite(l) ? -l * scale : std::numeric_limits<double>::infinity();
    };
    double best = std::numeric_limits<double>::infinity();
    for (auto& s : starts) {
      detail::require(s.size() == k, "a start with one value per parameter");
      calx::optimisation::NelderMeadOptions no;
      no.initial_step = {0.2};
      no.x_tolerance = 1e-10, no.f_tolerance = 1e-14;
      const auto nm = calx::optimisation::nelder_mead(f, s, no);
      std::vector<double> t = nm.x;
      double v = nm.value;
      bool converged = nm.converged;
      calx::optimisation::BfgsOptions bo;
      bo.gradient_tolerance = o.gradient_tolerance > 0 ? std::max(o.gradient_tolerance, 1e-9) : 1e-9;
      bo.max_iterations = o.max_iterations;
      const auto bf = calx::optimisation::bfgs(f, nm.x, bo);
      // BFGS's polish is kept unless it is worse beyond rounding: near the optimum the values tie to their last digits.
      if (std::isfinite(bf.value) && bf.value <= v + 1e-13 * (1 + std::fabs(v))) t = bf.x, v = std::min(v, bf.value), converged = bf.converged || converged;
      if (v < best) best = v, r.theta = t, r.converged = converged;
    }
    // Newton's steps on the gradient and Hessian by differences, while the gradient falls: values near the optimum tie to their
    // rounding, so a search on values stops about √ε from it, while the gradient by central differences is accurate to ε^{2/3}.
    auto gnorm = [&](const std::vector<double>& gr) {
      double m = 0;
      for (double a : gr) m = std::max(m, std::fabs(a) * scale);
      return m;
    };
    auto gr = calx::optimisation::gradient_by_differences(loglik, std::span<const double>(r.theta));
    Matrix H = calx::optimisation::hessian_by_differences(loglik, std::span<const double>(r.theta));
    for (int step = 0; step < 8; step++) {
      Matrix info(k);
      for (std::size_t a = 0; a < k; a++)
        for (std::size_t b = 0; b < k; b++) info(a, b) = -H(a, b);
      const auto c = calx::linalg::cholesky(info);
      if (!c.positive_definite) break;
      const auto dx = calx::linalg::cholesky_solve(c, gr);
      std::vector<double> next = r.theta;
      for (std::size_t a = 0; a < k; a++) next[a] += dx[a];
      if (!std::isfinite(loglik(next))) break;
      auto gn = calx::optimisation::gradient_by_differences(loglik, std::span<const double>(next));
      if (gnorm(gn) >= gnorm(gr)) break;
      r.theta = std::move(next), gr = std::move(gn);
      H = calx::optimisation::hessian_by_differences(loglik, std::span<const double>(r.theta));
    }
    r.converged = r.converged || gnorm(gr) <= 1e-8;
    r.log_likelihood = loglik(r.theta);
    r.information_positive = detail::information_inverse(H, r.covariance);
    if (!r.information_positive) r.covariance = Matrix(k);
    r.se = detail::standard_errors(r.covariance, r.information_positive, k);
  }
  r.aic = -2 * r.log_likelihood + 2 * static_cast<double>(k);
  if (q >= 1) r.sigma = std::exp(r.theta[0]);
  if (q == 2) r.shape = std::exp(r.theta[1]);
  auto eff = [&](std::span<const double> t) { return detail::effective(model, survey.transect, t, w); };
  r.effective = eff(r.theta);
  if (k > 0 && r.information_positive) {
    const auto grad = calx::optimisation::gradient_by_differences(eff, std::span<const double>(r.theta));
    double v = 0;
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t b = 0; b < k; b++) v += grad[a] * grad[b] * r.covariance(a, b);
    r.effective_se = std::sqrt(std::max(0.0, v));
  } else {
    r.effective_se = k == 0 ? 0 : not_a_number;
  }
  r.effective_cv = r.effective_se / r.effective;
  r.f0 = survey.transect == Transect::line ? 1 / r.effective : 2 * std::numbers::pi / r.effective;
  if (std::isfinite(w)) r.probability = survey.transect == Transect::line ? r.effective / w : r.effective / (std::numbers::pi * w * w);
  return r;
}

struct Density {
  double density = 0, se = 0, cv = 0, lower = 0, upper = 0, df = 0;
  double encounter_rate = 0, encounter_rate_se = 0;  // n/E
  double cv_encounter = 0, cv_detection = 0;
  std::size_t n = 0, transects = 0;
  double effort = 0;                           // E: total length, or total visits
  double covered_area = not_a_number;          // 2wL or (Σ visits) πw²
  double covered_abundance = not_a_number;     // n/P̂_a, the Horvitz–Thompson estimate in the covered area
};

inline Density density(const DistanceSurvey& survey, const DetectionFit& fit, double confidence = 0.95) {
  const std::size_t K = survey.effort.size();
  detail::require(K >= 1 && survey.on.size() == survey.distance.size(), "one transect for every detection, and the transects' efforts");
  detail::require(fit.truncation == survey.truncation && fit.transect == survey.transect, "the fit made on this survey's truncation and kind");
  std::vector<double> nk(K, 0.0);
  for (std::size_t i = 0; i < survey.distance.size(); i++) {
    detail::require(survey.on[i] < K, "transects numbered 0 … K − 1");
    if (survey.distance[i] <= survey.truncation) nk[survey.on[i]]++;
  }
  Density d;
  d.transects = K;
  for (std::size_t k = 0; k < K; k++) {
    detail::require(survey.effort[k] > 0, "positive efforts");
    d.effort += survey.effort[k];
    d.n += static_cast<std::size_t>(nk[k]);
  }
  const double n = static_cast<double>(d.n), E = d.effort, Kd = static_cast<double>(K);
  d.encounter_rate = n / E;
  if (K >= 2) {
    double s = 0;
    for (std::size_t k = 0; k < K; k++) {
      const double dev = nk[k] / survey.effort[k] - n / E;
      s += survey.effort[k] * survey.effort[k] * dev * dev;
    }
    d.encounter_rate_se = std::sqrt(Kd / (E * E * (Kd - 1)) * s);
  } else {
    d.encounter_rate_se = not_a_number;
  }
  d.cv_encounter = d.encounter_rate_se / d.encounter_rate;
  d.cv_detection = fit.effective_cv;
  d.density = survey.transect == Transect::line ? d.encounter_rate / (2 * fit.effective) : d.encounter_rate / fit.effective;
  d.cv = std::sqrt(d.cv_encounter * d.cv_encounter + d.cv_detection * d.cv_detection);
  d.se = d.cv * d.density;
  const double q = static_cast<double>(fit.parameters);
  const double denom = std::pow(d.cv_encounter, 4) / (Kd - 1) + (fit.parameters > 0 ? std::pow(d.cv_detection, 4) / (n - q) : 0.0);
  d.df = denom > 0 ? std::pow(d.cv, 4) / denom : not_a_number;
  const double t = std::isfinite(d.df) && d.df > 0 ? tyche::probability::StudentT{d.df}.upper_quantile((1 - confidence) / 2) : detail::z_two_sided(confidence);
  const double C = std::exp(t * std::sqrt(std::log(1 + d.cv * d.cv)));
  d.lower = d.density / C, d.upper = d.density * C;
  if (std::isfinite(survey.truncation)) {
    const double w = survey.truncation;
    d.covered_area = survey.transect == Transect::line ? 2 * w * E : E * std::numbers::pi * w * w;
    d.covered_abundance = n / fit.probability;
  }
  return d;
}

struct DistanceAnalysis {
  std::vector<DetectionFit> fits;  // one per model, in the order given
  std::vector<double> delta_aic;   // AIC − the least AIC
  std::size_t chosen = 0;          // the model of least AIC
  Density density;                 // with the chosen model
};

inline DistanceAnalysis distance_sampling(const DistanceSurvey& survey, std::span<const DetectionModel> models, double confidence = 0.95) {
  detail::require(!models.empty(), "one detection model at least");
  DistanceAnalysis a;
  for (const auto& m : models) a.fits.push_back(fit_detection(survey, m));
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < a.fits.size(); i++)
    if (a.fits[i].aic < best) best = a.fits[i].aic, a.chosen = i;
  for (const auto& f : a.fits) a.delta_aic.push_back(f.aic - best);
  a.density = density(survey, a.fits[a.chosen], confidence);
  return a;
}

}  // namespace peira::detection
