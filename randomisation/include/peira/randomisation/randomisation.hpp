// Inference for randomised experiments, from the design alone (docs/theory.md, §3).
//
// Each unit i has two potential outcomes Y_i(1), Y_i(0), fixed; the design is the only randomness. The estimand is the
// finite-population average effect τ = (1/n) Σ (Y_i(1) − Y_i(0)).
//
// - Neyman (1923): under complete randomisation the difference in means τ̂ is unbiased, with variance
//     Var(τ̂) = S₁²/n₁ + S₀²/n₀ − S_τ²/n,
//   S_τ² the variance of the unit effects, which no experiment sees. The estimator V̂ = s₁²/n₁ + s₀²/n₀ has expectation
//   S₁²/n₁ + S₀²/n₀ ≥ Var(τ̂), with equality exactly when the effect is constant (Imbens and Rubin 2015, Thm 6.2).
//   Blocks: the blocks' estimates weighted by their shares n_b/n, their variances by the squares of the shares.
// - Fisher (1935): the sharp null Y_i(1) = Y_i(0) + τ₀ for every unit fills in every potential outcome, so the law of any
//   statistic over the design is known. The statistic here is the difference in means of Y − τ₀ Z; its randomisation
//   distribution is enumerated when it has at most `max_exact` assignments, else drawn (p = (b + 1)/(B + 1), Phipson and
//   Smyth 2010, through tyche's testing).
//
// Limits: the intervals of Neyman are the normal's (large n₁, n₀); Fisher's test is exact for its sharp null, which is
// stronger than a null of no average effect.

#pragma once

#include <peira/design/design.hpp>

#include <tyche/probability/distributions.hpp>
#include <tyche/testing/permutation.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace peira::randomisation {

struct Neyman {
  double estimate = 0, variance = 0, se = 0, ci_lower = 0, ci_upper = 0;
  std::size_t n1 = 0, n0 = 0;
};

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::randomisation: ") + what);
}
struct Arm { double n = 0, mean = 0, var = 0; };
inline Arm arm(std::span<const double> y, std::span<const int> z, int which) {
  Arm a;
  for (std::size_t i = 0; i < y.size(); i++) if (z[i] == which) { a.n++; a.mean += y[i]; }
  a.mean /= a.n;
  for (std::size_t i = 0; i < y.size(); i++) if (z[i] == which) a.var += (y[i] - a.mean) * (y[i] - a.mean);
  a.var /= a.n - 1;
  return a;
}
inline double z_two_sided(double confidence) { return tyche::probability::Normal{0, 1}.upper_quantile((1 - confidence) / 2); }
}  // namespace detail

// Neyman's difference in means with its conservative variance (each arm at least two units).
inline Neyman difference_in_means(std::span<const double> y, std::span<const int> z, double confidence = 0.95) {
  detail::require(y.size() == z.size(), "one assignment per outcome");
  const auto a1 = detail::arm(y, z, 1), a0 = detail::arm(y, z, 0);
  detail::require(a1.n >= 2 && a0.n >= 2, "at least two units in each arm");
  Neyman r;
  r.n1 = static_cast<std::size_t>(a1.n), r.n0 = static_cast<std::size_t>(a0.n);
  r.estimate = a1.mean - a0.mean;
  r.variance = a1.var / a1.n + a0.var / a0.n;
  r.se = std::sqrt(r.variance);
  const double q = detail::z_two_sided(confidence);
  r.ci_lower = r.estimate - q * r.se, r.ci_upper = r.estimate + q * r.se;
  return r;
}

// Blocked (stratified) experiments: Σ_b (n_b/n) τ̂_b, with variance Σ_b (n_b/n)² V̂_b.
inline Neyman stratified(std::span<const double> y, std::span<const int> z, std::span<const int> block, double confidence = 0.95) {
  detail::require(y.size() == z.size() && y.size() == block.size(), "one assignment and one block per outcome");
  std::map<int, std::vector<std::size_t>> members;
  for (std::size_t i = 0; i < block.size(); i++) members[block[i]].push_back(i);
  Neyman r;
  const double n = static_cast<double>(y.size());
  for (const auto& [b, units] : members) {
    std::vector<double> yb;
    std::vector<int> zb;
    for (std::size_t i : units) { yb.push_back(y[i]); zb.push_back(z[i]); }
    const Neyman nb = difference_in_means(yb, zb, confidence);
    const double share = static_cast<double>(units.size()) / n;
    r.estimate += share * nb.estimate;
    r.variance += share * share * nb.variance;
    r.n1 += nb.n1, r.n0 += nb.n0;
  }
  r.se = std::sqrt(r.variance);
  const double q = detail::z_two_sided(confidence);
  r.ci_lower = r.estimate - q * r.se, r.ci_upper = r.estimate + q * r.se;
  return r;
}

enum class Alternative { two_sided, greater, less };

struct FisherTest {
  double observed = 0;      // the difference in means of Y − τ₀ Z
  double p_value = 1;
  std::size_t assignments = 0;  // enumerated, or drawn
  bool exact = false;
};

// Fisher's randomisation test of the sharp null of a constant effect τ₀, under complete randomisation of the observed n₁.
inline FisherTest fisher(std::span<const double> y, std::span<const int> z, double tau0 = 0, Alternative alternative = Alternative::two_sided,
                         std::size_t max_exact = 200000, std::size_t draws = 20000, std::uint64_t seed = 1) {
  detail::require(y.size() == z.size(), "one assignment per outcome");
  const std::size_t n = y.size();
  std::vector<double> y0(n);  // the outcomes under control, filled in by the sharp null
  std::size_t n1 = 0;
  for (std::size_t i = 0; i < n; i++) { y0[i] = y[i] - tau0 * z[i]; n1 += z[i] == 1; }
  detail::require(n1 >= 1 && n1 < n, "both arms non-empty");
  const auto statistic = [&](std::span<const int> w) {
    double s1 = 0, s0 = 0;
    for (std::size_t i = 0; i < n; i++) (w[i] ? s1 : s0) += y0[i];
    return s1 / static_cast<double>(n1) - s0 / static_cast<double>(n - n1);
  };
  FisherTest r;
  r.observed = statistic(z);
  // A small tolerance, so that assignments equal to the observed one in exact arithmetic count as extreme.
  const double tol = 1e-12 * (1 + std::fabs(r.observed));
  const auto extreme = [&](double t) {
    switch (alternative) {
      case Alternative::greater: return t >= r.observed - tol;
      case Alternative::less: return t <= r.observed + tol;
      default: return std::fabs(t) >= std::fabs(r.observed) - tol;
    }
  };
  if (peira::design::binomial_coefficient(n, n1) <= static_cast<double>(max_exact)) {
    std::size_t total = 0, hits = 0;
    peira::design::each_complete(n, n1, [&](const peira::design::Assignment& w) {
      total++;
      hits += extreme(statistic(w));
    });
    r.exact = true, r.assignments = total;
    r.p_value = static_cast<double>(hits) / static_cast<double>(total);
  } else {
    peira::design::Stream st(seed, 0x70e1u, 0x0001u, 0);
    std::size_t hits = 0;
    for (std::size_t b = 0; b < draws; b++) hits += extreme(statistic(peira::design::complete(n, n1, st)));
    r.assignments = draws;
    r.p_value = tyche::testing::permutation_p_value(hits, draws);
  }
  return r;
}

}  // namespace peira::randomisation
