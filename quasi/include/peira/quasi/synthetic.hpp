// Synthetic control: one treated unit compared with a weighted average of untreated units chosen to reproduce its past (Abadie and
// Gardeazabal 2003; Abadie, Diamond and Hainmueller 2010, J. Am. Stat. Assoc. 105, 493–505; Abadie 2021, J. Econ. Lit. 59, 391–425).
//
// THE WEIGHTS. With Y₁ the treated unit's outcomes and Y₀ the J donors', over the T₀ periods before treatment, the weights w lie in
// the simplex (w_j ≥ 0, Σ w_j = 1) and minimise ‖Y₁ − Y₀ w‖² over those periods: the outcomes are the predictors, each weighted
// alike (V = I, a declared choice; covariates and a chosen V are not taken yet). The problem is a convex quadratic over the simplex,
// solved by accelerated projected gradient (FISTA; Beck and Teboulle 2009) with the exact Euclidean projection onto the simplex
// (Duchi, Shalev-Shwartz, Singer and Chandra 2008), to a stated tolerance on the change of w. Being in the simplex, the synthetic
// unit interpolates and never extrapolates; its fit before treatment says whether it can stand for the treated unit at all.
//
// THE EFFECT. The gaps Y₁t − Σ_j w_j Y_jt after T₀; the root mean squared prediction errors before and after, and their ratio.
//
// INFERENCE IN SPACE (ADH 2010, §2.4; Abadie et al. 2015). Each donor in turn is treated as if treated, with the other donors as its
// pool, and its ratio computed; the p-value is the rank of the treated unit's ratio among the J + 1, from the largest: the share of
// units whose post-to-pre ratio is at least as large. With J donors the smallest attainable p-value is 1/(J + 1).

#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

#include "common.hpp"

namespace peira::quasi {

// The Euclidean projection of v onto the simplex {w ≥ 0, Σ w = 1} (Duchi et al. 2008, fig. 1).
inline std::vector<double> project_simplex(const std::vector<double>& v) {
  std::vector<double> u = v;
  std::sort(u.begin(), u.end(), std::greater<>());
  double sum = 0, theta = 0;
  for (std::size_t j = 0; j < u.size(); j++) {
    sum += u[j];
    const double t = (sum - 1) / static_cast<double>(j + 1);
    if (u[j] - t > 0) theta = t;
  }
  std::vector<double> w(v.size());
  for (std::size_t j = 0; j < v.size(); j++) w[j] = std::max(0.0, v[j] - theta);
  return w;
}

struct SyntheticControl {
  std::vector<double> weights;          // per donor
  std::vector<double> synthetic, gap;   // per period
  double rmspe_pre = 0, rmspe_post = 0, ratio = 0;
  std::size_t iterations = 0;
  // Inference in space.
  std::vector<double> placebo_ratio;    // per donor treated as if treated
  double p_value = std::numeric_limits<double>::quiet_NaN();
};

namespace detail {
// Weights of a target series over donors' series (rows: donors) on the first t0 periods.
inline std::vector<double> simplex_weights(const std::vector<double>& target, const std::vector<std::vector<double>>& donors, std::size_t t0, double tolerance,
                                           std::size_t max_iterations, std::size_t& iterations) {
  const std::size_t J = donors.size();
  // A = Y₀ᵀY₀ and b = Y₀ᵀY₁ on the pre-periods; the gradient of ‖Y₁ − Y₀w‖² is 2(Aw − b).
  std::vector<double> A(J * J, 0.0), b(J, 0.0);
  for (std::size_t i = 0; i < J; i++) {
    for (std::size_t t = 0; t < t0; t++) b[i] += donors[i][t] * target[t];
    for (std::size_t j = 0; j < J; j++)
      for (std::size_t t = 0; t < t0; t++) A[i * J + j] += donors[i][t] * donors[j][t];
  }
  // The Lipschitz constant 2 λ_max(A), bounded by power iteration (and the trace as a fallback).
  double lambda = 0;
  {
    std::vector<double> v(J, 1.0 / std::sqrt(static_cast<double>(J))), Av(J);
    for (int k = 0; k < 200; k++) {
      double norm = 0;
      for (std::size_t i = 0; i < J; i++) {
        Av[i] = 0;
        for (std::size_t j = 0; j < J; j++) Av[i] += A[i * J + j] * v[j];
        norm += Av[i] * Av[i];
      }
      norm = std::sqrt(norm);
      if (!(norm > 0)) break;
      lambda = norm;
      for (std::size_t i = 0; i < J; i++) v[i] = Av[i] / norm;
    }
    double trace = 0;
    for (std::size_t i = 0; i < J; i++) trace += A[i * J + i];
    lambda = std::min(trace, 1.05 * lambda + 1e-300);
    if (!(lambda > 0)) lambda = 1;
  }
  const double step = 1 / (2 * lambda);
  std::vector<double> w(J, 1.0 / static_cast<double>(J)), z = w, previous = w, g(J);
  double t = 1;
  iterations = 0;
  for (; iterations < max_iterations; iterations++) {
    for (std::size_t i = 0; i < J; i++) {
      g[i] = -b[i];
      for (std::size_t j = 0; j < J; j++) g[i] += A[i * J + j] * z[j];
      g[i] *= 2;
    }
    std::vector<double> y(J);
    for (std::size_t i = 0; i < J; i++) y[i] = z[i] - step * g[i];
    previous = w;
    w = project_simplex(y);
    const double tn = (1 + std::sqrt(1 + 4 * t * t)) / 2;
    double change = 0;
    for (std::size_t i = 0; i < J; i++) {
      z[i] = w[i] + (t - 1) / tn * (w[i] - previous[i]);
      change = std::max(change, std::fabs(w[i] - previous[i]));
    }
    t = tn;
    if (change < tolerance && iterations > 10) break;
  }
  return w;
}

inline void fill(SyntheticControl& s, const std::vector<double>& target, const std::vector<std::vector<double>>& donors, std::size_t t0) {
  const std::size_t T = target.size();
  s.synthetic.assign(T, 0.0), s.gap.assign(T, 0.0);
  double pre = 0, post = 0;
  for (std::size_t t = 0; t < T; t++) {
    for (std::size_t j = 0; j < donors.size(); j++) s.synthetic[t] += s.weights[j] * donors[j][t];
    s.gap[t] = target[t] - s.synthetic[t];
    (t < t0 ? pre : post) += s.gap[t] * s.gap[t];
  }
  s.rmspe_pre = std::sqrt(pre / static_cast<double>(t0));
  s.rmspe_post = T > t0 ? std::sqrt(post / static_cast<double>(T - t0)) : 0;
  s.ratio = s.rmspe_pre > 0 ? s.rmspe_post / s.rmspe_pre : std::numeric_limits<double>::infinity();
}
}  // namespace detail

// The synthetic control of `treated` (its outcomes over T periods) from `donors` (J series of T), treated from period t0 on;
// `placebos` runs the inference in space.
inline SyntheticControl synthetic_control(const std::vector<double>& treated, const std::vector<std::vector<double>>& donors, std::size_t t0, bool placebos = true,
                                          double tolerance = 1e-12, std::size_t max_iterations = 200000) {
  detail::require(donors.size() >= 2, "two donors at least");
  detail::require(t0 >= 1 && t0 <= treated.size(), "a pre-treatment period at least");
  for (const auto& d : donors) detail::require(d.size() == treated.size(), "every donor observed in every period");
  SyntheticControl s;
  s.weights = detail::simplex_weights(treated, donors, t0, tolerance, max_iterations, s.iterations);
  detail::fill(s, treated, donors, t0);
  if (placebos) {
    std::size_t at_least = 1;  // the treated unit itself
    for (std::size_t j = 0; j < donors.size(); j++) {
      std::vector<std::vector<double>> pool;
      for (std::size_t k = 0; k < donors.size(); k++)
        if (k != j) pool.push_back(donors[k]);
      SyntheticControl p;
      std::size_t it = 0;
      p.weights = detail::simplex_weights(donors[j], pool, t0, tolerance, max_iterations, it);
      detail::fill(p, donors[j], pool, t0);
      s.placebo_ratio.push_back(p.ratio);
      at_least += p.ratio >= s.ratio;
    }
    s.p_value = static_cast<double>(at_least) / static_cast<double>(donors.size() + 1);
  }
  return s;
}

}  // namespace peira::quasi
