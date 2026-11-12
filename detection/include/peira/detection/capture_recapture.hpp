// Closed populations by capture and recapture (docs/theory.md, §10.4).
//
// FULL DATA AND MECHANISM. A population of N individuals, unknown, does not change while it is sampled (closure). On each of t
// occasions individuals are caught, marked if new and released; individual h on occasion j is caught with probability p_j. The
// full data are the N capture histories; the mechanism hides every history of zeros, so what is seen is the M histories with at
// least one capture. N = M + (the number never caught), and the marks identify p from the recaptures.
//
// IDENTIFICATION STATED. (i) Closure: no births, deaths or movement between occasions (untestable from the histories alone in
// these models). (ii) Marks are kept and read without error. (iii) Within an occasion every individual has the same probability
// of capture, and captures are independent (untestable: individuals that are easier to catch make N̂ too small; models with
// heterogeneity are not written here). (iv) Model Mt lets p vary by occasion; M0 holds it constant.
//
// - lincoln_petersen: n1 marked on the first occasion, n2 caught on the second, m2 of them marked: N̂ = n1 n2/m2 (Petersen 1896;
//   Lincoln 1930), with the variance n1 n2 (n1 − m2)(n2 − m2)/m2³, the delta method's under the hypergeometric law of m2 given
//   n1 and n2. It needs m2 > 0 and is biased upward in small samples.
// - chapman: N̂ = (n1 + 1)(n2 + 1)/(m2 + 1) − 1 (Chapman 1951), whose bias under the hypergeometric law is exactly
//   −(n2 + 1) C(N − n1, n2 + 1)/C(N, n2), zero when n1 + n2 ≥ N (theory note, Theorem 10.5; checked by enumeration), with Seber's
//   variance (n1 + 1)(n2 + 1)(n1 − m2)(n2 − m2)/((m2 + 1)²(m2 + 2)) (Seber 1970; Seber
//   1982, The Estimation of Animal Abundance, 2nd ed., §3.1.1).
//   Both give the normal interval and Chao's (1987, Biometrics 43, 783–791) log-normal interval for the number never caught,
//   f̂0 = N̂ − M: [M + f̂0/C, M + f̂0 C], C = exp(z √log(1 + v/f̂0²)), which never falls below the M seen.
// - schnabel: t occasions with C_t caught, R_t of them marked before, M_t the marked at large before occasion t:
//   N̂ = Σ C_t M_t / Σ R_t (Schnabel 1938), and Chapman's form Σ C_t M_t / (Σ R_t + 1); the interval treats Σ R_t as Poisson,
//   with its exact limits from the χ² quantiles (Seber 1982). With two occasions it is Lincoln–Petersen.
// - closed_m0, closed_mt: maximum likelihood over N (Otis, Burnham, White and Anderson 1978, Statistical inference from capture
//   data on closed animal populations, Wildlife Monographs 62; Darroch 1958, Biometrika 45, 343–359, for Mt). The likelihood is
//   L(N, p) ∝ N!/(N − M)! Π_j p_j^{n_j} (1 − p_j)^{N − n_j}, n_j the captures on occasion j (M0: p_j = p), up to the factor
//   1/Π_ω x_ω! of the histories' counts, the same for every model of the same histories, so that AICs compare. For fixed N the
//   probabilities are p̂_j = n_j/N (M0: p̂ = Σ n_j/(tN)), and N̂ maximises the profile ℓ(N, p̂(N)) over the integers N ≥ M,
//   found by stepping N upward; the interval is the profile likelihood's, {N : 2(ℓ(N̂) − ℓ(N)) ≤ χ²_1(confidence)}, which
//   never reaches below the M seen and is not symmetric about N̂ as a normal interval is. With no recapture the profile
//   increases without bound: N̂ and the upper limit are then +∞ (reported when the search reaches max_abundance, by default
//   1000 M + 1000).
//   Darroch's equation 1 − M/N = Π_j (1 − n_j/N) is the stationarity of this profile with log(N!/(N − M)!)′ replaced by
//   −log(1 − M/N); it says N = M/(1 − Π_j(1 − p̂_j)), the Horvitz–Thompson form (§10.1 of the theory note).

#pragma once

#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include "common.hpp"

namespace peira::detection {

struct Petersen {
  double abundance = 0, se = 0;
  double lower = 0, upper = 0;                // Chao's log-normal interval
  double normal_lower = 0, normal_upper = 0;  // N̂ ± z se
  double seen = 0;                            // M = n1 + n2 − m2
};

namespace detail {
inline Petersen petersen_intervals(double N, double v, double M, double confidence) {
  Petersen r;
  r.abundance = N, r.se = std::sqrt(v), r.seen = M;
  const double z = z_two_sided(confidence), f0 = N - M;
  r.normal_lower = N - z * r.se, r.normal_upper = N + z * r.se;
  if (f0 > 0) {
    const double C = std::exp(z * std::sqrt(std::log(1 + v / (f0 * f0))));
    r.lower = M + f0 / C, r.upper = M + f0 * C;
  } else {
    r.lower = r.upper = N;
  }
  return r;
}
inline void require_petersen(double n1, double n2, double m2) {
  require(n1 >= 0 && n2 >= 0 && m2 >= 0 && m2 <= std::min(n1, n2), "0 ≤ m2 ≤ min(n1, n2)");
}
}  // namespace detail

inline Petersen lincoln_petersen(double n1, double n2, double m2, double confidence = 0.95) {
  detail::require_petersen(n1, n2, m2);
  detail::require(m2 > 0, "at least one recapture");
  return detail::petersen_intervals(n1 * n2 / m2, n1 * n2 * (n1 - m2) * (n2 - m2) / (m2 * m2 * m2), n1 + n2 - m2, confidence);
}

inline Petersen chapman(double n1, double n2, double m2, double confidence = 0.95) {
  detail::require_petersen(n1, n2, m2);
  const double N = (n1 + 1) * (n2 + 1) / (m2 + 1) - 1;
  const double v = (n1 + 1) * (n2 + 1) * (n1 - m2) * (n2 - m2) / ((m2 + 1) * (m2 + 1) * (m2 + 2));
  return detail::petersen_intervals(N, v, n1 + n2 - m2, confidence);
}

// Capture histories: individuals × occasions, row-major, 1 caught and 0 not. Rows of zeros (never caught) are not allowed.
struct CaptureHistories {
  std::size_t individuals = 0, occasions = 0;
  std::vector<int> y;

  CaptureHistories() = default;
  CaptureHistories(std::size_t n, std::size_t t) : individuals(n), occasions(t), y(n * t, 0) {}
  int& operator()(std::size_t h, std::size_t j) { return y[h * occasions + j]; }
  int operator()(std::size_t h, std::size_t j) const { return y[h * occasions + j]; }
};

struct CaptureSummary {
  std::vector<double> caught;         // n_j (C_j): caught on occasion j
  std::vector<double> recaptured;     // R_j: of them, marked before
  std::vector<double> marked_before;  // M_j: marked at large before occasion j (M_1 = 0)
  double distinct = 0;                // M_{t+1}: individuals caught at least once
  double captures = 0;                // Σ n_j
};

inline CaptureSummary summarise(const CaptureHistories& h) {
  detail::require(h.y.size() == h.individuals * h.occasions && h.occasions >= 1, "a table of individuals by occasions");
  CaptureSummary s;
  s.caught.assign(h.occasions, 0.0), s.recaptured.assign(h.occasions, 0.0), s.marked_before.assign(h.occasions, 0.0);
  std::vector<int> first(h.individuals, -1);
  for (std::size_t i = 0; i < h.individuals; i++) {
    for (std::size_t j = 0; j < h.occasions; j++) {
      const int v = h(i, j);
      detail::require(v == 0 || v == 1, "captures 0 or 1");
      if (v && first[i] < 0) first[i] = static_cast<int>(j);
    }
    detail::require(first[i] >= 0, "every individual caught at least once");
  }
  for (std::size_t i = 0; i < h.individuals; i++)
    for (std::size_t j = 0; j < h.occasions; j++)
      if (h(i, j)) s.caught[j]++, s.recaptured[j] += static_cast<int>(j) > first[i];
  double marked = 0;
  for (std::size_t j = 0; j < h.occasions; j++) {
    s.marked_before[j] = marked;
    marked += s.caught[j] - s.recaptured[j];
    s.captures += s.caught[j];
  }
  s.distinct = marked;
  return s;
}

struct Schnabel {
  double abundance = 0;           // Σ C M / Σ R
  double abundance_chapman = 0;   // Σ C M / (Σ R + 1)
  double lower = 0, upper = 0;    // from Σ R as Poisson
  double recaptures = 0, sum_caught_marked = 0;
};

// caught C_t and recaptured R_t on each occasion; the marked at large M_t follow (M_1 = 0, M_{t+1} = M_t + C_t − R_t).
inline Schnabel schnabel(std::span<const double> caught, std::span<const double> recaptured, double confidence = 0.95) {
  detail::require(caught.size() == recaptured.size() && caught.size() >= 2, "two occasions at least, with their recaptures");
  Schnabel r;
  double marked = 0;
  for (std::size_t t = 0; t < caught.size(); t++) {
    detail::require(recaptured[t] >= 0 && recaptured[t] <= caught[t] && recaptured[t] <= marked, "0 ≤ R_t ≤ min(C_t, M_t)");
    r.sum_caught_marked += caught[t] * marked;
    r.recaptures += recaptured[t];
    marked += caught[t] - recaptured[t];
  }
  const double R = r.recaptures, a = 1 - confidence;
  r.abundance = R > 0 ? r.sum_caught_marked / R : std::numeric_limits<double>::infinity();
  r.abundance_chapman = r.sum_caught_marked / (R + 1);
  // The exact Poisson limits for the mean of Σ R: χ²_{2R}(a/2)/2 and χ²_{2R+2}(1 − a/2)/2.
  const double mean_lo = R > 0 ? tyche::probability::ChiSquared{2 * R}.quantile(a / 2) / 2 : 0;
  const double mean_hi = tyche::probability::ChiSquared{2 * R + 2}.upper_quantile(a / 2) / 2;
  r.lower = r.sum_caught_marked / mean_hi;
  r.upper = mean_lo > 0 ? r.sum_caught_marked / mean_lo : std::numeric_limits<double>::infinity();
  return r;
}

// log L(N, p) for model Mt (p one per occasion; M0 is p constant), up to the histories' constant: log N!/(N − M)! +
// Σ_j [n_j log p_j + (N − n_j) log(1 − p_j)], for real N ≥ M.
inline double closed_log_likelihood(double N, std::span<const double> p, std::span<const double> caught, double distinct) {
  detail::require(p.size() == caught.size(), "one probability per occasion");
  if (N < distinct) return -std::numeric_limits<double>::infinity();
  double l = calx::special_functions::log_gamma(N + 1) - calx::special_functions::log_gamma(N - distinct + 1);
  for (std::size_t j = 0; j < p.size(); j++) {
    if (caught[j] > 0) l += caught[j] * std::log(p[j]);
    if (N - caught[j] > 0) l += (N - caught[j]) * std::log1p(-p[j]);
  }
  return l;
}

struct ClosedPopulation {
  double abundance = 0;          // N̂, an integer ≥ M (or +∞)
  double lower = 0, upper = 0;   // the profile likelihood interval, integers (upper +∞ when open within max_abundance)
  std::vector<double> p;         // at N̂: one (M0) or one per occasion (Mt)
  double probability_seen = 0;   // 1 − Π_j (1 − p̂_j): the probability of being caught at least once
  double log_likelihood = 0, aic = 0;
  std::size_t parameters = 0;
  double distinct = 0;           // M
};

namespace detail {
// The profile over the integers N ≥ M of a function ℓ(N), stepping upward until past the maximum and past the interval.
template <class Profile>
inline ClosedPopulation closed_profile(Profile&& profile, double M, double confidence, double max_abundance) {
  ClosedPopulation r;
  r.distinct = M;
  const double cut = tyche::probability::ChiSquared{1}.upper_quantile(1 - confidence) / 2;
  const double top = max_abundance > 0 ? max_abundance : 1000 * M + 1000;
  std::vector<double> values;  // ℓ(M + k)
  double best = -std::numeric_limits<double>::infinity(), argbest = M;
  for (double N = M; N <= top; N++) {
    const double l = profile(N);
    values.push_back(l);
    if (l > best) best = l, argbest = N;
    if (N > argbest && best - l > cut + 1) break;
  }
  const bool open = best - values.back() <= cut;
  r.abundance = argbest;
  r.log_likelihood = best;
  r.lower = M + static_cast<double>(values.size());
  r.upper = M;
  for (std::size_t k = 0; k < values.size(); k++)
    if (best - values[k] <= cut) r.lower = std::min(r.lower, M + static_cast<double>(k)), r.upper = std::max(r.upper, M + static_cast<double>(k));
  if (open) {
    r.upper = std::numeric_limits<double>::infinity();
    if (argbest >= top) r.abundance = std::numeric_limits<double>::infinity();
  }
  return r;
}
}  // namespace detail

// Model M0 from its sufficient statistics: Σ n_j captures in all, M distinct, t occasions.
inline ClosedPopulation closed_m0(double captures, double distinct, std::size_t occasions, double confidence = 0.95, double max_abundance = 0) {
  detail::require(occasions >= 2 && distinct >= 1 && captures >= distinct && captures <= distinct * static_cast<double>(occasions),
                  "two occasions at least, and M ≤ Σ n_j ≤ tM");
  const double t = static_cast<double>(occasions);
  auto profile = [&](double N) {
    const double p = captures / (t * N);
    double l = calx::special_functions::log_gamma(N + 1) - calx::special_functions::log_gamma(N - distinct + 1) + captures * std::log(p);
    if (t * N - captures > 0) l += (t * N - captures) * std::log1p(-p);
    return l;
  };
  auto r = detail::closed_profile(profile, distinct, confidence, max_abundance);
  const double p = std::isfinite(r.abundance) ? captures / (t * r.abundance) : 0;
  r.p = {p};
  r.probability_seen = 1 - std::pow(1 - p, t);
  r.parameters = 2;
  r.aic = -2 * r.log_likelihood + 4;
  return r;
}

// Model Mt from its sufficient statistics: n_j caught on each occasion, M distinct.
inline ClosedPopulation closed_mt(std::span<const double> caught, double distinct, double confidence = 0.95, double max_abundance = 0) {
  detail::require(caught.size() >= 2 && distinct >= 1, "two occasions at least, and one individual");
  double sum = 0;
  for (double n : caught) {
    detail::require(n >= 0 && n <= distinct, "0 ≤ n_j ≤ M");
    sum += n;
  }
  detail::require(sum >= distinct, "Σ n_j ≥ M");
  std::vector<double> p(caught.size());
  auto profile = [&](double N) {
    for (std::size_t j = 0; j < caught.size(); j++) p[j] = caught[j] / N;
    return closed_log_likelihood(N, p, caught, distinct);
  };
  auto r = detail::closed_profile(profile, distinct, confidence, max_abundance);
  r.p.resize(caught.size());
  double none = 1;
  for (std::size_t j = 0; j < caught.size(); j++) {
    r.p[j] = std::isfinite(r.abundance) ? caught[j] / r.abundance : 0;
    none *= 1 - r.p[j];
  }
  r.probability_seen = 1 - none;
  r.parameters = caught.size() + 1;
  r.aic = -2 * r.log_likelihood + 2 * static_cast<double>(r.parameters);
  return r;
}

inline ClosedPopulation closed_m0(const CaptureHistories& h, double confidence = 0.95, double max_abundance = 0) {
  const auto s = summarise(h);
  return closed_m0(s.captures, s.distinct, h.occasions, confidence, max_abundance);
}
inline ClosedPopulation closed_mt(const CaptureHistories& h, double confidence = 0.95, double max_abundance = 0) {
  const auto s = summarise(h);
  return closed_mt(s.caught, s.distinct, confidence, max_abundance);
}

}  // namespace peira::detection
