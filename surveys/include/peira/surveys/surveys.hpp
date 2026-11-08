// Samples drawn by a design of strata and clusters, and the estimates they allow (docs/theory.md, §4).
//
// A sample's units carry their weights w_i = 1/π_i, their stratum h and their primary sampling unit (PSU); a stratum's
// PSUs are drawn independently of the other strata's, at a sampling fraction f_h (0: with replacement, or a fraction
// too small to matter). Every estimate is a weighted total, or a smooth function of totals, and its variance is the
// ultimate-cluster estimator
//     V̂ = Σ_h (1 − f_h) n_h/(n_h − 1) Σ_{a ∈ h} (z_ha − z̄_h)²,
// z_ha the weighted total over PSU a of the estimate's linearised variable (Wolter 2007, Introduction to Variance
// Estimation, 2nd ed., §1.4 and ch. 6; Särndal, Swensson and Wretman 1992, §5.5–5.7), with n_h − 1 = PSUs − strata
// degrees of freedom. It is exactly unbiased for stratified simple random sampling of units without replacement, and for
// the Hansen–Hurwitz estimator of PSUs drawn with replacement, and approximately so otherwise.
//
// - total, mean, ratio: Σ w y; the ratio R̂ = Σ w y / Σ w x linearised as u = (y − R̂ x)/Σ w x; the mean the ratio
//   to x = 1.
// - stratified_srswor: the design of a stratified simple random sample of units, from the strata's population sizes.
// - rake: weights adjusted by iterative proportional fitting until their sums over each margin's categories meet known
//   totals (Deming and Stephan 1940); post-stratification is raking on one margin, reached in one step.
//
// Limits: estimates for domains (subpopulations) are made by setting y = 0 outside them, not by subsetting the sample;
// a stratum with one PSU has no variance estimate and is refused.

#pragma once

#include <tyche/probability/distributions.hpp>

#include <cmath>
#include <cstddef>
#include <map>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace peira::surveys {

struct Design {
  std::vector<double> weight;  // w_i = 1/π_i
  std::vector<int> stratum;    // h; one stratum for all if empty
  std::vector<int> psu;        // the unit's primary sampling unit within its stratum; each unit its own if empty
  std::vector<double> fraction;  // f_h of the unit's stratum (the PSUs' sampling fraction); 0 for all if empty
};

struct Estimate {
  double value = 0, se = 0, df = 0, ci_lower = 0, ci_upper = 0;
};

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::surveys: ") + what);
}

// The ultimate-cluster variance of Σ w u, and its degrees of freedom (PSUs − strata).
inline std::pair<double, double> variance(std::span<const double> u, const Design& d) {
  const std::size_t n = u.size();
  require(d.weight.size() == n, "one weight per unit");
  require(d.stratum.empty() || d.stratum.size() == n, "one stratum per unit");
  require(d.psu.empty() || d.psu.size() == n, "one PSU per unit");
  require(d.fraction.empty() || d.fraction.size() == n, "one sampling fraction per unit");
  std::map<int, std::map<long long, double>> totals;  // stratum → PSU → Σ w u
  std::map<int, double> fraction;
  for (std::size_t i = 0; i < n; i++) {
    const int h = d.stratum.empty() ? 0 : d.stratum[i];
    const long long a = d.psu.empty() ? static_cast<long long>(i) : d.psu[i];
    totals[h][a] += d.weight[i] * u[i];
    fraction[h] = d.fraction.empty() ? 0 : d.fraction[i];
  }
  double v = 0, psus = 0;
  for (const auto& [h, by_psu] : totals) {
    const double nh = static_cast<double>(by_psu.size());
    require(nh >= 2, "at least two PSUs in every stratum");
    double mean = 0;
    for (const auto& [a, t] : by_psu) mean += t / nh;
    double ss = 0;
    for (const auto& [a, t] : by_psu) ss += (t - mean) * (t - mean);
    v += (1 - fraction[h]) * nh / (nh - 1) * ss;
    psus += nh;
  }
  return {v, psus - static_cast<double>(totals.size())};
}

inline Estimate finish(double value, std::pair<double, double> vdf, double confidence) {
  Estimate e;
  e.value = value, e.se = std::sqrt(vdf.first), e.df = vdf.second;
  const double q = tyche::probability::StudentT{e.df}.upper_quantile((1 - confidence) / 2);
  e.ci_lower = value - q * e.se, e.ci_upper = value + q * e.se;
  return e;
}
}  // namespace detail

inline Estimate total(std::span<const double> y, const Design& d, double confidence = 0.95) {
  double t = 0;
  for (std::size_t i = 0; i < y.size(); i++) t += d.weight[i] * y[i];
  return detail::finish(t, detail::variance(y, d), confidence);
}

inline Estimate ratio(std::span<const double> y, std::span<const double> x, const Design& d, double confidence = 0.95) {
  detail::require(y.size() == x.size(), "one x per y");
  double ty = 0, tx = 0;
  for (std::size_t i = 0; i < y.size(); i++) { ty += d.weight[i] * y[i]; tx += d.weight[i] * x[i]; }
  detail::require(tx != 0, "a non-zero total of x");
  const double R = ty / tx;
  std::vector<double> u(y.size());
  for (std::size_t i = 0; i < y.size(); i++) u[i] = (y[i] - R * x[i]) / tx;
  return detail::finish(R, detail::variance(u, d), confidence);
}

inline Estimate mean(std::span<const double> y, const Design& d, double confidence = 0.95) {
  const std::vector<double> one(y.size(), 1.0);
  return ratio(y, one, d, confidence);
}

// A stratified simple random sample of units without replacement: population sizes N_h of the strata, the sample's
// strata; weights N_h/n_h and fractions n_h/N_h.
inline Design stratified_srswor(std::span<const int> stratum, const std::map<int, double>& population) {
  std::map<int, double> n;
  for (int h : stratum) n[h]++;
  Design d;
  d.stratum.assign(stratum.begin(), stratum.end());
  for (int h : stratum) {
    const auto it = population.find(h);
    detail::require(it != population.end() && it->second >= n[h], "a population size N_h ≥ n_h for every stratum");
    d.weight.push_back(it->second / n[h]);
    d.fraction.push_back(n[h] / it->second);
  }
  return d;
}

struct Raking {
  std::vector<double> weight;
  int iterations = 0;
  bool converged = false;
  double gap = 0;  // the largest |Σ w over a category − its total| at the end
};

// Iterative proportional fitting of the weights to the totals of each margin: margins[m][i] is unit i's category in
// margin m, totals[m] the known total of each category.
inline Raking rake(std::span<const double> weight, const std::vector<std::vector<int>>& margins, const std::vector<std::map<int, double>>& totals,
                   int max_iterations = 200, double tolerance = 1e-10) {
  detail::require(margins.size() == totals.size(), "a total for every margin");
  Raking r;
  r.weight.assign(weight.begin(), weight.end());
  const std::size_t n = r.weight.size();
  for (const auto& m : margins) detail::require(m.size() == n, "a category in every margin for every unit");
  for (r.iterations = 1; r.iterations <= max_iterations; r.iterations++) {
    for (std::size_t m = 0; m < margins.size(); m++) {
      std::map<int, double> sum;
      for (std::size_t i = 0; i < n; i++) sum[margins[m][i]] += r.weight[i];
      for (std::size_t i = 0; i < n; i++) {
        const auto it = totals[m].find(margins[m][i]);
        detail::require(it != totals[m].end(), "a total for every category met in the sample");
        r.weight[i] *= it->second / sum[margins[m][i]];
      }
    }
    r.gap = 0;
    for (std::size_t m = 0; m < margins.size(); m++) {
      std::map<int, double> sum;
      for (std::size_t i = 0; i < n; i++) sum[margins[m][i]] += r.weight[i];
      for (const auto& [c, t] : totals[m]) r.gap = std::max(r.gap, std::fabs(sum[c] - t) / std::max(1.0, std::fabs(t)));
    }
    if (r.gap <= tolerance) { r.converged = true; break; }
  }
  return r;
}

}  // namespace peira::surveys
