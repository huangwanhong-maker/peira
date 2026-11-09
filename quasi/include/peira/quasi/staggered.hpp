// Treatment adopted at different times: group-time average effects and their aggregations (Callaway and Sant'Anna 2021,
// "Difference-in-differences with multiple time periods", J. Econometrics 225, 200–230), and the conventional event study by two-way
// fixed effects beside them.
//
// THE DESIGN. Units i are observed in periods t; each is first treated in period G_i and stays treated (an absorbing treatment), or
// never (G_i = never). The units first treated in g form group g.
//
// GROUP-TIME EFFECTS. Under parallel trends with a comparison group C (units never treated, or also those not yet treated by t),
//   ATT(g, t) = E[Y_t − Y_{g−1} | G = g] − E[Y_t − Y_{g−1} | C],
// the average effect at t on the units first treated in g, against the universal base period g − 1. It is identified for t ≥ g;
// for t < g it is a placebo, zero under parallel trends (the pre-trends of an event study). Each mean is over the units observed
// in both t and g − 1. No covariates, so parallel trends are unconditional (CS 2021, §2–3).
//
// AGGREGATIONS (CS 2021, §3.1). With N_g the size of group g:
// - simple: the average of ATT(g, t) over t ≥ g, weighted by N_g;
// - by exposure e = t − g (the event study): θ(e) = Σ_g N_g ATT(g, g + e) / Σ_g N_g, over the groups in which g + e is observed.
// Their standard errors are from a nonparametric bootstrap over units (the whole series of a unit drawn together), from a calx
// stream given by the caller; CS 2021 use a multiplier bootstrap, asymptotically equivalent.
//
// UNITS REACHED BEFORE THEY ARE TREATED. With interference a unit may be touched by others' treatments (a spillover) before, or
// without, its own. An optional `reached` gives each unit's first period touched by any treatment, no later than its first treated
// period; a unit then serves as a comparison only while not yet reached (never reached, under never_treated). Without it, reached is
// the first treated period, and the estimator is the one above (interference.hpp, peira's decision 0006).
//
// TWO-WAY FIXED EFFECTS. The static regression of Y on unit and period effects and the treatment indicator D_it = [t ≥ G_i], and
// the dynamic one on relative-period indicators (e = −1 left out), estimated by tyche's within estimator. With effects that vary
// across groups or over exposure, their coefficients are weighted sums of the ATT(g, t) with weights that may be negative
// (Goodman-Bacon 2021; Sun and Abraham 2021; de Chaisemartin and D'Haultfœuille 2020): reported here for comparison.

#pragma once

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <utility>
#include <vector>

#include "common.hpp"

namespace peira::quasi {

constexpr std::int64_t never = std::numeric_limits<std::int64_t>::max();

enum class Comparison { never_treated, not_yet_treated };

struct GroupTime {
  std::int64_t group = 0, period = 0;
  double att = 0, se = std::numeric_limits<double>::quiet_NaN();
  std::size_t treated = 0, comparison = 0;  // units in each mean
};

struct StaggeredFit {
  std::vector<GroupTime> effects;                       // every (g, t) estimable
  double simple = 0, simple_se = std::numeric_limits<double>::quiet_NaN();
  std::map<std::int64_t, double> by_exposure, by_exposure_se;  // θ(e) and its standard error
  std::map<std::int64_t, std::size_t> group_size;      // N_g
};

namespace detail {

// y[unit][period] for a long panel; units and periods by their labels.
struct Wide {
  std::vector<std::int64_t> units, periods;
  std::vector<std::map<std::int64_t, double>> y;        // per unit: period → outcome
  std::vector<std::int64_t> first;                      // per unit: G_i
  std::vector<std::int64_t> reached;                    // per unit: the first period touched by any treatment, ≤ G_i
};

inline Wide wide(std::span<const double> y, std::span<const std::int64_t> unit, std::span<const std::int64_t> period, std::span<const std::int64_t> first,
                 std::span<const std::int64_t> reached = {}) {
  require(unit.size() == y.size() && period.size() == y.size() && first.size() == y.size(), "a unit, a period and a first treated period for every observation");
  require(reached.empty() || reached.size() == y.size(), "a first period reached for every observation, or none");
  Wide w;
  std::map<std::int64_t, std::size_t> index;
  std::set<std::int64_t> periods;
  for (std::size_t i = 0; i < y.size(); i++) {
    const std::int64_t ri = reached.empty() ? first[i] : reached[i];
    require(ri <= first[i], "a unit reached no later than it is treated");
    auto [it, fresh] = index.emplace(unit[i], w.units.size());
    if (fresh) w.units.push_back(unit[i]), w.y.emplace_back(), w.first.push_back(first[i]), w.reached.push_back(ri);
    require(w.first[it->second] == first[i] && w.reached[it->second] == ri, "one first treated and one first reached period per unit");
    w.y[it->second][period[i]] = y[i];
    periods.insert(period[i]);
  }
  w.periods.assign(periods.begin(), periods.end());
  return w;
}

// ATT(g, t) on the units drawn (with multiplicity), with its counts.
inline GroupTime group_time(const Wide& w, const std::vector<std::size_t>& draw, std::int64_t g, std::int64_t t, std::int64_t base, Comparison c) {
  GroupTime r;
  r.group = g, r.period = t;
  double st = 0, sc = 0;
  for (const std::size_t i : draw) {
    const auto& yi = w.y[i];
    const auto a = yi.find(t), b = yi.find(base);
    if (a == yi.end() || b == yi.end()) continue;
    const std::int64_t gi = w.first[i], ri = w.reached[i];
    if (gi == g) st += a->second - b->second, r.treated++;
    else if (ri == never || (c == Comparison::not_yet_treated && ri > std::max(t, base))) sc += a->second - b->second, r.comparison++;
  }
  r.att = (r.treated > 0 && r.comparison > 0) ? st / static_cast<double>(r.treated) - sc / static_cast<double>(r.comparison) : std::numeric_limits<double>::quiet_NaN();
  return r;
}

struct Estimates {
  std::vector<GroupTime> effects;
  double simple = 0;
  std::map<std::int64_t, double> by_exposure;
};

inline Estimates estimate(const Wide& w, const std::vector<std::size_t>& draw, Comparison c, std::map<std::int64_t, std::size_t>* sizes = nullptr) {
  std::map<std::int64_t, std::size_t> size;
  for (const std::size_t i : draw)
    if (w.first[i] != never) size[w.first[i]]++;
  if (sizes) *sizes = size;
  Estimates e;
  double num = 0, den = 0;
  std::map<std::int64_t, std::pair<double, double>> expo;
  for (const auto& [g, n] : size) {
    const std::int64_t base = g - 1;
    for (const std::int64_t t : w.periods) {
      if (t == base) continue;
      const auto r = group_time(w, draw, g, t, base, c);
      if (!std::isfinite(r.att)) continue;
      e.effects.push_back(r);
      if (t >= g) num += static_cast<double>(n) * r.att, den += static_cast<double>(n);
      auto& x = expo[t - g];
      x.first += static_cast<double>(n) * r.att, x.second += static_cast<double>(n);
    }
  }
  e.simple = den > 0 ? num / den : std::numeric_limits<double>::quiet_NaN();
  for (const auto& [k, v] : expo) e.by_exposure[k] = v.first / v.second;
  return e;
}

}  // namespace detail

// Group-time effects and their aggregations; `first` is each observation's unit's first treated period (never for none), and
// `reached` (optional) its first period touched by any treatment. Bootstrap standard errors from `replications` draws of units (0 for
// none) on the stream given.
inline StaggeredFit staggered(std::span<const double> y, std::span<const std::int64_t> unit, std::span<const std::int64_t> period, std::span<const std::int64_t> first,
                              Comparison comparison = Comparison::never_treated, std::size_t replications = 0, calx::random::Stream* stream = nullptr,
                              std::span<const std::int64_t> reached = {}) {
  const auto w = detail::wide(y, unit, period, first, reached);
  std::vector<std::size_t> all(w.units.size());
  for (std::size_t i = 0; i < all.size(); i++) all[i] = i;
  StaggeredFit fit;
  const auto e = detail::estimate(w, all, comparison, &fit.group_size);
  fit.effects = e.effects;
  fit.simple = e.simple;
  fit.by_exposure = e.by_exposure;
  if (replications > 0) {
    detail::require(stream != nullptr, "a stream for the bootstrap");
    std::map<std::pair<std::int64_t, std::int64_t>, std::pair<double, double>> gt;
    std::map<std::int64_t, std::pair<double, double>> ex;
    std::pair<double, double> simple{0, 0};
    std::map<std::pair<std::int64_t, std::int64_t>, double> gt_n;
    std::map<std::int64_t, double> ex_n;
    double simple_n = 0;
    std::vector<std::size_t> draw(all.size());
    for (std::size_t b = 0; b < replications; b++) {
      for (auto& d : draw) d = static_cast<std::size_t>(stream->below(all.size()));
      const auto r = detail::estimate(w, draw, comparison);
      for (const auto& x : r.effects) {
        auto& s = gt[{x.group, x.period}];
        s.first += x.att, s.second += x.att * x.att, gt_n[{x.group, x.period}] += 1;
      }
      for (const auto& [k, v] : r.by_exposure)
        if (std::isfinite(v)) ex[k].first += v, ex[k].second += v * v, ex_n[k] += 1;
      if (std::isfinite(r.simple)) simple.first += r.simple, simple.second += r.simple * r.simple, simple_n += 1;
    }
    const auto sd = [](const std::pair<double, double>& s, double n) { return n > 1 ? std::sqrt(std::max(0.0, (s.second - s.first * s.first / n) / (n - 1))) : std::numeric_limits<double>::quiet_NaN(); };
    for (auto& x : fit.effects) x.se = sd(gt[{x.group, x.period}], gt_n[{x.group, x.period}]);
    for (const auto& [k, v] : fit.by_exposure) fit.by_exposure_se[k] = sd(ex[k], ex_n[k]);
    fit.simple_se = sd(simple, simple_n);
  }
  return fit;
}

// The static two-way fixed-effects regression on D_it = [t ≥ G_i], errors clustered by unit.
inline WithinFit staggered_two_way(std::span<const double> y, std::span<const std::int64_t> unit, std::span<const std::int64_t> period, std::span<const std::int64_t> first) {
  Rectangular X(y.size(), 1);
  for (std::size_t i = 0; i < y.size(); i++) X(i, 0) = first[i] != never && period[i] >= first[i] ? 1.0 : 0.0;
  PanelOptions o;
  o.effects = PanelEffects::two_way;
  o.covariance = PanelCovariance::cluster;
  return within(X, y, unit, period, o);
}

// The dynamic two-way regression on relative-period indicators e ∈ [lo, hi], e ≠ −1 (periods beyond are binned at the ends),
// never-treated units having none; its coefficients by e.
struct EventStudy {
  std::map<std::int64_t, double> coefficient, se;
  WithinFit fit;
};
inline EventStudy event_study_two_way(std::span<const double> y, std::span<const std::int64_t> unit, std::span<const std::int64_t> period, std::span<const std::int64_t> first,
                                      std::int64_t lo, std::int64_t hi) {
  detail::require(lo < -1 && hi >= 0, "a window with leads before −1 and lags from 0");
  std::vector<std::int64_t> es;
  for (std::int64_t e = lo; e <= hi; e++)
    if (e != -1) es.push_back(e);
  Rectangular X(y.size(), es.size());
  for (std::size_t i = 0; i < y.size(); i++) {
    if (first[i] == never) continue;
    const std::int64_t e = std::clamp(period[i] - first[i], lo, hi);
    for (std::size_t k = 0; k < es.size(); k++) X(i, k) = es[k] == e ? 1.0 : 0.0;
  }
  PanelOptions o;
  o.effects = PanelEffects::two_way;
  o.covariance = PanelCovariance::cluster;
  EventStudy s;
  s.fit = within(X, y, unit, period, o);
  for (std::size_t k = 0; k < es.size(); k++) s.coefficient[es[k]] = s.fit.fit.beta[k], s.se[es[k]] = s.fit.fit.se[k];
  return s;
}

}  // namespace peira::quasi
