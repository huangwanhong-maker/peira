// The power of a design and the sample it needs, by the normal approximation to the law of its estimator.
//
// - Means: two arms of n₁ and n₀ units, a difference δ, a common standard deviation σ, a two-sided test at level α.
//   The difference in means has standard error σ √(1/n₁ + 1/n₀), so the power is
//     Φ(|δ|/se − z_{1−α/2}) + Φ(−|δ|/se − z_{1−α/2}),
//   and equal arms need n = 2 σ² (z_{1−α/2} + z_{1−β})² / δ² each (Cohen 1988, §2.4; the t law's correction adds about
//   one unit at small n, and is not applied).
// - Proportions: p₁ against p₀, without continuity correction (Fleiss, Levin and Paik 2003, Statistical Methods for Rates
//   and Proportions, 3rd ed., eq. 4.14): n = [z_{1−α/2} √(2 p̄ q̄) + z_{1−β} √(p₁q₁ + p₀q₀)]² / (p₁ − p₀)².
// - Clusters of m units with intraclass correlation ρ multiply the variance by the design effect 1 + (m − 1) ρ
//   (Kish 1965, Survey Sampling, §5.4), and so the units needed.
//
// Limits: the approximation is the normal's; for small samples a design's power should be checked by simulating it.

#pragma once

#include <tyche/probability/distributions.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace peira::design {

namespace detail {
inline double z_upper(double q) { return tyche::probability::Normal{0, 1}.upper_quantile(q); }
inline double phi(double x) { return tyche::probability::Normal{0, 1}.cdf(x); }
}  // namespace detail

inline double power_means(double delta, double sigma, double n1, double n0, double alpha = 0.05) {
  if (!(sigma > 0 && n1 > 0 && n0 > 0 && alpha > 0 && alpha < 1)) throw std::invalid_argument("peira::design::power_means: σ, n₁, n₀ > 0 and α in (0, 1)");
  const double se = sigma * std::sqrt(1 / n1 + 1 / n0), z = detail::z_upper(alpha / 2), d = std::fabs(delta) / se;
  return detail::phi(d - z) + detail::phi(-d - z);
}

inline std::size_t n_per_arm_means(double delta, double sigma, double power = 0.8, double alpha = 0.05) {
  if (!(delta != 0 && sigma > 0 && power > 0 && power < 1)) throw std::invalid_argument("peira::design::n_per_arm_means: δ ≠ 0, σ > 0, power in (0, 1)");
  const double z = detail::z_upper(alpha / 2) + detail::z_upper(1 - power);
  return static_cast<std::size_t>(std::ceil(2 * sigma * sigma * z * z / (delta * delta) - 1e-12));
}

inline double power_proportions(double p1, double p0, double n_per_arm, double alpha = 0.05) {
  if (!(p1 > 0 && p1 < 1 && p0 > 0 && p0 < 1 && p1 != p0 && n_per_arm > 0)) throw std::invalid_argument("peira::design::power_proportions: p₁ ≠ p₀ in (0, 1)");
  const double pbar = (p1 + p0) / 2, z = detail::z_upper(alpha / 2);
  const double s0 = std::sqrt(2 * pbar * (1 - pbar)), s1 = std::sqrt(p1 * (1 - p1) + p0 * (1 - p0));
  return detail::phi((std::fabs(p1 - p0) * std::sqrt(n_per_arm) - z * s0) / s1);
}

inline std::size_t n_per_arm_proportions(double p1, double p0, double power = 0.8, double alpha = 0.05) {
  if (!(p1 > 0 && p1 < 1 && p0 > 0 && p0 < 1 && p1 != p0)) throw std::invalid_argument("peira::design::n_per_arm_proportions: p₁ ≠ p₀ in (0, 1)");
  const double pbar = (p1 + p0) / 2;
  const double a = detail::z_upper(alpha / 2) * std::sqrt(2 * pbar * (1 - pbar)) + detail::z_upper(1 - power) * std::sqrt(p1 * (1 - p1) + p0 * (1 - p0));
  return static_cast<std::size_t>(std::ceil(a * a / ((p1 - p0) * (p1 - p0)) - 1e-12));
}

inline double design_effect(double m, double rho) {
  if (!(m >= 1)) throw std::invalid_argument("peira::design::design_effect: clusters of at least one unit");
  return 1 + (m - 1) * rho;
}

}  // namespace peira::design
