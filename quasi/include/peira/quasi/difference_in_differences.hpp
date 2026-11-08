// Differences in differences.
//
// THE 2 × 2 DESIGN. Units are treated or not, observed before and after the treatment; with ȳ_gp
// the mean outcome of group g ∈ {control, treated} in period p ∈ {pre, post},
//   δ̂ = (ȳ_treated,post − ȳ_treated,pre) − (ȳ_control,post − ȳ_control,pre),
// the average effect of the treatment on the treated under parallel trends (Card and Krueger 1994;
// Angrist and Pischke 2009, ch. 5). It is the coefficient of treated × post in the OLS of y on an
// intercept, treated, post and their product, whose standard error is computed by tyche's OLS:
// clustered by unit when unit labels are given (Bertrand, Duflo and Mullainathan 2004), HC1
// otherwise. δ̂ equals the difference of the four means exactly, the regression being saturated.
//
// TWO-WAY FIXED EFFECTS. With many units and periods, y_it = α_i + γ_t + δ D_it + x_it β + ε_it,
// D_it the treatment indicator (or intensity), estimated by the two-way within estimator
// (panel.hpp) with standard errors clustered by unit. In a balanced 2 × 2 panel it reproduces δ̂
// above (checked).
//
// Caveat. With treatment adopted at different times and effects that vary over time, the two-way
// regression averages comparisons with negative weights and need not estimate any average effect
// (Goodman-Bacon 2021; de Chaisemartin and D'Haultfœuille 2020): staggered.hpp holds the
// estimators robust to this (Callaway and Sant'Anna 2021) and event studies. Moved from nomos in
// peira's stage 2 (peira's decision 0004).

#pragma once


#include <cstdint>
#include <span>
#include <vector>

#include "common.hpp"

namespace peira::quasi {

struct DifferenceInDifferences {
  double estimate = 0;                  // δ̂
  double se = 0, t = 0, p = 0, ci_lower = 0, ci_upper = 0;
  double mean[2][2] = {{0, 0}, {0, 0}};  // ȳ[treated][post]
  double count[2][2] = {{0, 0}, {0, 0}};
  tr::LinearFit regression;             // y on 1, treated, post, treated × post
};

// The 2 × 2 estimator; treated and post are 0 or 1 for each observation; units (optional) cluster the errors.
inline DifferenceInDifferences difference_in_differences(std::span<const double> y, std::span<const int> treated, std::span<const int> post,
                                                         std::span<const std::int64_t> units = {}, double confidence = 0.95) {
  const std::size_t n = y.size();
  detail::require(treated.size() == n && post.size() == n && (units.empty() || units.size() == n), "a group and a period for every observation");
  DifferenceInDifferences d;
  Rectangular X(n, 4);
  for (std::size_t i = 0; i < n; i++) {
    detail::require((treated[i] == 0 || treated[i] == 1) && (post[i] == 0 || post[i] == 1), "indicators of 0 or 1");
    d.mean[treated[i]][post[i]] += y[i];
    d.count[treated[i]][post[i]] += 1;
    X(i, 0) = 1;
    X(i, 1) = treated[i];
    X(i, 2) = post[i];
    X(i, 3) = treated[i] * post[i];
  }
  for (int g = 0; g < 2; g++)
    for (int p = 0; p < 2; p++) {
      detail::require(d.count[g][p] > 0, "observations in each of the four cells");
      d.mean[g][p] /= d.count[g][p];
    }
  tr::OlsOptions o;
  o.confidence = confidence;
  if (units.empty()) {
    o.covariance = tr::Covariance::hc1;
  } else {
    o.covariance = tr::Covariance::cluster;
    o.clusters.assign(units.begin(), units.end());
  }
  d.regression = tr::ols(X, y, o);
  d.estimate = d.regression.beta[3];
  d.se = d.regression.se[3];
  d.t = d.regression.t[3];
  d.p = d.regression.p[3];
  d.ci_lower = d.regression.ci_lower[3];
  d.ci_upper = d.regression.ci_upper[3];
  return d;
}

// The two-way fixed-effects regression of y on the treatment D (its first coefficient) and the
// controls X (optional), with unit and period effects and errors clustered by unit.
inline WithinFit two_way_fixed_effects(std::span<const double> y, std::span<const double> treatment, std::span<const std::int64_t> unit,
                                       std::span<const std::int64_t> period, const Rectangular& controls = {}, double confidence = 0.95) {
  const std::size_t n = y.size();
  detail::require(treatment.size() == n && (controls.cols == 0 || controls.rows == n), "a treatment and controls for every observation");
  Rectangular X(n, 1 + controls.cols);
  for (std::size_t i = 0; i < n; i++) {
    X(i, 0) = treatment[i];
    for (std::size_t j = 0; j < controls.cols; j++) X(i, j + 1) = controls(i, j);
  }
  PanelOptions o;
  o.effects = PanelEffects::two_way;
  o.covariance = PanelCovariance::cluster;
  o.confidence = confidence;
  return within(X, y, unit, period, o);
}

}  // namespace peira::quasi
