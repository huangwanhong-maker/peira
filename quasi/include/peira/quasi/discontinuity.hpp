// Regression discontinuity: a treatment assigned (sharply) or encouraged (fuzzily) where a running variable crosses a cutoff
// (Thistlethwaite and Campbell 1960; Hahn, Todd and van der Klaauw 2001; Imbens and Lemieux 2008, J. Econometrics 142, 615–635).
//
// SHARP. D = [x ≥ c]. If the potential outcomes' means are continuous at c, the effect at the cutoff is the jump
//   τ = lim_{x↓c} E[Y | x] − lim_{x↑c} E[Y | x],
// estimated by local linear regression on each side (Fan and Gijbels 1996; Hahn et al. 2001): the weighted least squares of y on
// 1, D, (x − c) and D (x − c) with kernel weights K((x − c)/h), τ̂ its coefficient on D. The kernel is triangular (the default,
// optimal at a boundary; Cheng, Fan and Marron 1997), Epanechnikov or uniform; h the bandwidth. Standard errors are the WLS
// sandwich, heteroskedasticity-robust (HC1), conditional on h: no correction for the bias of local linear fits (Calonico,
// Cattaneo and Titiunik 2014), so a bandwidth chosen to minimise the mean squared error leaves intervals that undercover.
//
// FUZZY. Treatment d jumps in probability at c. Under monotonicity and continuity, τ_F = (jump in E[Y | x]) / (jump in E[d | x]),
// the local average effect on compliers at the cutoff (Hahn et al. 2001). It is the instrumental-variables estimate with D as
// the instrument for d and 1, (x − c), D (x − c) exogenous, on the √K-weighted data (instrumental.hpp): just identified, it is the
// ratio of the two local linear jumps exactly (checked).
//
// BANDWIDTH. Given, or chosen by leave-one-out cross-validation for the boundary (Ludwig and Miller 2007; Imbens and Lemieux
// 2008, §4.5): each observation within the share `window` of the closest to c on its side is predicted from the observations on its
// own side between it and h away from the cutoff (the side a boundary estimate sees), by a local linear fit with the kernel; the h
// of a grid with the least mean squared error is chosen.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <vector>

#include "common.hpp"
#include "instrumental.hpp"

namespace peira::quasi {

enum class Kernel { triangular, epanechnikov, uniform };

inline double kernel_weight(Kernel k, double u) {
  const double a = std::fabs(u);
  if (a >= 1) return 0;
  switch (k) {
    case Kernel::triangular: return 1 - a;
    case Kernel::epanechnikov: return 0.75 * (1 - a * a);
    case Kernel::uniform: return 0.5;
  }
  return 0;
}

struct Discontinuity {
  double estimate = 0, se = 0, t = 0, p = 0, ci_lower = 0, ci_upper = 0;
  double bandwidth = 0;
  std::size_t left = 0, right = 0;  // observations with positive weight on each side
  double jump_outcome = 0, jump_treatment = 1;  // fuzzy: the two jumps; sharp: τ̂ and 1
};

namespace detail {
struct Local {
  Rectangular X;
  std::vector<double> y, d, w;
  std::size_t left = 0, right = 0;
};
inline Local local(std::span<const double> x, std::span<const double> y, std::span<const double> d, double c, double h, Kernel k) {
  require(x.size() == y.size() && (d.empty() || d.size() == y.size()), "a running value and an outcome (and a treatment) for every observation");
  require(h > 0, "a bandwidth h > 0");
  Local L;
  std::vector<std::size_t> keep;
  for (std::size_t i = 0; i < x.size(); i++) {
    const double w = kernel_weight(k, (x[i] - c) / h);
    if (w > 0) keep.push_back(i), L.w.push_back(w), (x[i] >= c ? L.right : L.left)++;
  }
  require(L.left >= 3 && L.right >= 3, "three observations at least within the bandwidth on each side");
  L.X = Rectangular(keep.size(), 4);
  for (std::size_t r = 0; r < keep.size(); r++) {
    const std::size_t i = keep[r];
    const double D = x[i] >= c ? 1.0 : 0.0, u = x[i] - c;
    L.X(r, 0) = 1, L.X(r, 1) = D, L.X(r, 2) = u, L.X(r, 3) = D * u;
    L.y.push_back(y[i]);
    if (!d.empty()) L.d.push_back(d[i]);
  }
  return L;
}
}  // namespace detail

// The sharp design's effect at the cutoff c with bandwidth h.
inline Discontinuity sharp_discontinuity(std::span<const double> x, std::span<const double> y, double c, double h, Kernel k = Kernel::triangular, double confidence = 0.95) {
  const auto L = detail::local(x, y, {}, c, h, k);
  tr::OlsOptions o;
  o.weights = L.w;
  o.covariance = tr::Covariance::hc1;
  o.confidence = confidence;
  const auto f = tr::ols(L.X, L.y, o);
  Discontinuity r;
  r.estimate = f.beta[1], r.se = f.se[1], r.t = f.t[1], r.p = f.p[1], r.ci_lower = f.ci_lower[1], r.ci_upper = f.ci_upper[1];
  r.bandwidth = h, r.left = L.left, r.right = L.right, r.jump_outcome = r.estimate;
  return r;
}

// The fuzzy design's effect on compliers at the cutoff: d instrumented by D on the √K-weighted local sample.
inline Discontinuity fuzzy_discontinuity(std::span<const double> x, std::span<const double> y, std::span<const double> d, double c, double h, Kernel k = Kernel::triangular,
                                         double confidence = 0.95) {
  detail::require(d.size() == y.size(), "a treatment for every observation");
  const auto L = detail::local(x, y, d, c, h, k);
  const std::size_t n = L.y.size();
  Rectangular W(n, 3), E(n, 1), Z(n, 1);
  std::vector<double> ys(n);
  for (std::size_t i = 0; i < n; i++) {
    const double s = std::sqrt(L.w[i]);
    W(i, 0) = s * L.X(i, 0), W(i, 1) = s * L.X(i, 2), W(i, 2) = s * L.X(i, 3);
    E(i, 0) = s * L.d[i];
    Z(i, 0) = s * L.X(i, 1);
    ys[i] = s * L.y[i];
  }
  IvOptions o;
  o.covariance = IvCovariance::hc1;
  o.confidence = confidence;
  const auto f = two_stage_least_squares(ys, W, E, Z, o);
  Discontinuity r;
  r.estimate = f.beta[3], r.se = f.se[3], r.t = f.t[3], r.p = f.p[3], r.ci_lower = f.ci_lower[3], r.ci_upper = f.ci_upper[3];
  r.bandwidth = h, r.left = L.left, r.right = L.right;
  tr::OlsOptions plain;
  plain.weights = L.w;
  r.jump_outcome = tr::ols(L.X, L.y, plain).beta[1];
  r.jump_treatment = tr::ols(L.X, L.d, plain).beta[1];
  return r;
}

// The bandwidth of `grid` least in leave-one-out error at the boundary (see the header); NaN errors where too few points.
struct Bandwidth {
  double chosen = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> grid, error;
};
inline Bandwidth cross_validated_bandwidth(std::span<const double> x, std::span<const double> y, double c, const std::vector<double>& grid, double window = 0.5,
                                           Kernel k = Kernel::triangular) {
  detail::require(x.size() == y.size() && !grid.empty() && window > 0 && window <= 1, "points, a grid of bandwidths and a window in (0, 1]");
  // The observations predicted: on each side, the share `window` closest to c.
  std::vector<double> lefts, rights;
  for (const double v : x) (v < c ? lefts : rights).push_back(v);
  std::sort(lefts.begin(), lefts.end()), std::sort(rights.begin(), rights.end());
  const double lo = lefts.empty() ? c : lefts[static_cast<std::size_t>((1 - window) * static_cast<double>(lefts.size() - 1))];
  const double hi = rights.empty() ? c : rights[static_cast<std::size_t>(window * static_cast<double>(rights.size() - 1))];
  Bandwidth b;
  b.grid = grid;
  double best = std::numeric_limits<double>::infinity();
  for (const double h : grid) {
    double sum = 0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < x.size(); i++) {
      const bool left = x[i] < c;
      if (left ? x[i] < lo : x[i] > hi) continue;
      // A one-sided local linear fit at x_i from the points on its side, farther from c than x_i, within h.
      double s0 = 0, s1 = 0, s2 = 0, t0 = 0, t1 = 0;
      std::size_t used = 0;
      for (std::size_t j = 0; j < x.size(); j++) {
        if (j == i || (x[j] < c) != left) continue;
        const double u = x[j] - x[i];
        if (left ? (u >= 0 || -u >= h) : (u <= 0 || u >= h)) continue;
        const double w = kernel_weight(k, u / h);
        if (w <= 0) continue;
        s0 += w, s1 += w * u, s2 += w * u * u, t0 += w * y[j], t1 += w * u * y[j], used++;
      }
      const double det = s0 * s2 - s1 * s1;
      if (used < 3 || !(std::fabs(det) > 1e-12 * std::max(1.0, s0 * s2))) continue;
      const double fit = (s2 * t0 - s1 * t1) / det;  // the intercept at u = 0
      sum += (y[i] - fit) * (y[i] - fit), count++;
    }
    const double e = count > 0 ? sum / static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
    b.error.push_back(e);
    if (std::isfinite(e) && e < best) best = e, b.chosen = h;
  }
  return b;
}

}  // namespace peira::quasi
