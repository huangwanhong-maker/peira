// Missing values, and inference from several imputations of them (docs/theory.md, §6).
//
// A response indicator R reveals Y where R = 1. Missing at random (MAR) means P(R | Y, X) = P(R | X): given what is
// always seen, whether a value is seen does not depend on the value. Under MAR, imputations drawn from the predictive law
// of the missing values given the seen ones, m times, and analysed each as if complete, combine by Rubin's rules into
// valid inference (Rubin 1987, Multiple Imputation for Nonresponse in Surveys, ch. 3).
//
// - pool: Q̄ = mean of the estimates, Ū = mean of their variances, B = their variance between imputations,
//   T = Ū + (1 + 1/m) B; the degrees of freedom of Barnard and Rubin (1999), which tend to Rubin's (m − 1)/γ² as the
//   complete data's degrees of freedom grow, γ = (1 + 1/m) B / T the fraction of information missing (van Buuren 2018,
//   Flexible Imputation of Missing Data, 2nd ed., §2.3).
// - approximate_bayesian_bootstrap: within each class of X, the seen values resampled with replacement to a donor pool
//   of their number, from which the missing values are drawn (Rubin and Schenker 1986); a proper imputation under MAR
//   given the class, the donor pool's resampling carrying the uncertainty of the class's law.
//
// Limits: MAR cannot be tested from the data; under missingness not at random the imputations are biased, and only a
// sensitivity analysis (planned) speaks to it.

#pragma once

#include <peira/design/design.hpp>

#include <tyche/probability/distributions.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace peira::missing {

struct Pooled {
  double estimate = 0, within = 0, between = 0, total = 0, se = 0, df = 0, fraction_missing = 0, ci_lower = 0, ci_upper = 0;
  std::size_t imputations = 0;
};

// m ≥ 2 completed-data estimates and their variances; complete_df the degrees of freedom of the complete-data analysis
// (infinity for a large-sample analysis).
inline Pooled pool(std::span<const double> estimates, std::span<const double> variances, double complete_df = std::numeric_limits<double>::infinity(),
                   double confidence = 0.95) {
  const std::size_t m = estimates.size();
  if (m < 2 || variances.size() != m) throw std::invalid_argument("peira::missing::pool: at least two imputations, each with its variance");
  const double M = static_cast<double>(m);
  Pooled p;
  p.imputations = m;
  for (std::size_t l = 0; l < m; l++) { p.estimate += estimates[l] / M; p.within += variances[l] / M; }
  for (std::size_t l = 0; l < m; l++) p.between += (estimates[l] - p.estimate) * (estimates[l] - p.estimate) / (M - 1);
  p.total = p.within + (1 + 1 / M) * p.between;
  p.se = std::sqrt(p.total);
  const double gamma = p.total > 0 ? (1 + 1 / M) * p.between / p.total : 0;
  p.fraction_missing = gamma;
  const double nu_m = gamma > 0 ? (M - 1) / (gamma * gamma) : std::numeric_limits<double>::infinity();
  if (std::isinf(complete_df)) p.df = nu_m;
  else {
    const double nu_obs = (complete_df + 1) / (complete_df + 3) * complete_df * (1 - gamma);
    p.df = 1 / (1 / nu_m + 1 / nu_obs);
  }
  const double q = std::isinf(p.df) ? tyche::probability::Normal{0, 1}.upper_quantile((1 - confidence) / 2)
                                    : tyche::probability::StudentT{p.df}.upper_quantile((1 - confidence) / 2);
  p.ci_lower = p.estimate - q * p.se, p.ci_upper = p.estimate + q * p.se;
  return p;
}

// y with its seen values (observed[i] = 1) kept and its missing ones drawn, class by class. A class with no seen value
// cannot be imputed and is refused.
inline std::vector<double> approximate_bayesian_bootstrap(std::span<const double> y, std::span<const int> observed, std::span<const int> klass,
                                                          peira::design::Stream& s) {
  const std::size_t n = y.size();
  if (observed.size() != n || klass.size() != n) throw std::invalid_argument("peira::missing::approximate_bayesian_bootstrap: one indicator and one class per value");
  std::map<int, std::vector<double>> seen;
  for (std::size_t i = 0; i < n; i++) if (observed[i]) seen[klass[i]].push_back(y[i]);
  std::map<int, std::vector<double>> donors;
  for (const auto& [c, values] : seen) {
    auto& d = donors[c];
    for (std::size_t k = 0; k < values.size(); k++) d.push_back(values[static_cast<std::size_t>(s.below(values.size()))]);
  }
  std::vector<double> out(y.begin(), y.end());
  for (std::size_t i = 0; i < n; i++) {
    if (observed[i]) continue;
    const auto it = donors.find(klass[i]);
    if (it == donors.end()) throw std::invalid_argument("peira::missing::approximate_bayesian_bootstrap: a class with no seen value");
    out[i] = it->second[static_cast<std::size_t>(s.below(it->second.size()))];
  }
  return out;
}

}  // namespace peira::missing
