// Cohorts seen in time: person-time and events on the Lexis diagram, rates with exact intervals, and the standardisation of rates
// (peira's stage 3, begun for demography; the infrastructure's ADR 0085).
//
// THE LEXIS DIAGRAM (Lexis 1875; Keiding 1990, J. R. Stat. Soc. A 153, 1–40). A life born at calendar time b and followed from entry
// to exit is the segment {(t, t − b) : entry ≤ t < exit} in the plane of calendar time and age, at 45°. A grid of age bands and
// period bands cuts it into pieces: each piece's length is that life's person-time in the cell, and an event at exit is counted in
// the cell holding the exit. Summed over lives, person-time T_ap and events D_ap give the occurrence–exposure rates D_ap / T_ap, the
// Poisson likelihood's maximum when the hazard is constant within cells (Breslow and Day 1987, ch. 2). Time outside the grid is
// left out and reported. Split exactly: the pieces' lengths sum to the follow-up within the grid (checked).
//
// RATES. A rate D/T with an exact interval (Garwood 1936): [χ²_{α/2}(2D)/(2T), χ²_{1−α/2}(2D + 2)/(2T)], the Poisson means whose tails
// at D are α/2 (checked by inverting the Poisson's distribution function), the lower bound 0 when D = 0.
//
// STANDARDISATION. Rates compared across populations of different structure (by age, say) are standardised:
// - direct: Σ_a w_a r_a / Σ_a w_a with a standard population's weights w_a, its variance Σ_a w_a² D_a / T_a² / (Σ w_a)² under Poisson
//   counts, and a normal interval (a declared approximation; Fay and Feuer 1997 give an exact-tailed one);
// - indirect: the standardised ratio SMR = D / E of observed to expected events, E = Σ_a T_a R_a with the standard's rates R_a, its
//   interval Garwood's for D divided by E (Breslow and Day 1987, ch. 2).

#pragma once

#include <tyche/probability/distributions.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace peira::cohorts {

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::cohorts: ") + what);
}
}  // namespace detail

// A life: born at `birth`, followed from `entry` to `exit` (calendar times, the same unit as the grid), the event at exit or not.
struct Life {
  double birth = 0, entry = 0, exit = 0;
  bool event = false;
};

// Age bands [age[k], age[k + 1]) and period bands [period[j], period[j + 1]); edges increasing.
struct Grid {
  std::vector<double> age, period;
};

struct Lexis {
  std::vector<std::vector<double>> person_time, events;  // [age band][period band]
  double outside = 0;                                     // person-time outside the grid
  double events_outside = 0;
};

// Each life's follow-up split along the grid's age and period edges.
inline Lexis lexis(const std::vector<Life>& lives, const Grid& g) {
  detail::require(g.age.size() >= 2 && g.period.size() >= 2, "two edges at least in each direction");
  detail::require(std::is_sorted(g.age.begin(), g.age.end()) && std::is_sorted(g.period.begin(), g.period.end()), "increasing edges");
  const std::size_t A = g.age.size() - 1, P = g.period.size() - 1;
  Lexis L;
  L.person_time.assign(A, std::vector<double>(P, 0.0));
  L.events.assign(A, std::vector<double>(P, 0.0));
  const auto band = [](const std::vector<double>& e, double v) -> std::ptrdiff_t {
    if (v < e.front() || v >= e.back()) return -1;
    return std::upper_bound(e.begin(), e.end(), v) - e.begin() - 1;
  };
  const auto band_before = [](const std::vector<double>& e, double v) -> std::ptrdiff_t {
    if (v <= e.front() || v > e.back()) return -1;
    return std::lower_bound(e.begin(), e.end(), v) - e.begin() - 1;
  };
  for (const Life& life : lives) {
    detail::require(life.exit >= life.entry && life.entry >= life.birth, "birth ≤ entry ≤ exit");
    // The times at which the segment crosses an age edge (t = b + a_k) or a period edge (t = p_j), between entry and exit.
    std::vector<double> cut{life.entry, life.exit};
    for (const double a : g.age)
      if (life.birth + a > life.entry && life.birth + a < life.exit) cut.push_back(life.birth + a);
    for (const double p : g.period)
      if (p > life.entry && p < life.exit) cut.push_back(p);
    std::sort(cut.begin(), cut.end());
    for (std::size_t k = 0; k + 1 < cut.size(); k++) {
      const double t0 = cut[k], t1 = cut[k + 1], mid = (t0 + t1) / 2;
      if (!(t1 > t0)) continue;
      const auto a = band(g.age, mid - life.birth), p = band(g.period, mid);
      if (a < 0 || p < 0) L.outside += t1 - t0;
      else L.person_time[static_cast<std::size_t>(a)][static_cast<std::size_t>(p)] += t1 - t0;
    }
    if (life.event) {
      // The event at exit, in the cell holding the instant just before it: bands (e_k, e_{k+1}] for this lookup.
      const auto a = band_before(g.age, life.exit - life.birth), p = band_before(g.period, life.exit);
      if (a < 0 || p < 0) L.events_outside += 1;
      else L.events[static_cast<std::size_t>(a)][static_cast<std::size_t>(p)] += 1;
    }
  }
  return L;
}

struct Rate {
  double events = 0, exposure = 0, rate = 0, lower = 0, upper = 0;
};

// D/T with Garwood's exact interval at the given confidence.
inline Rate rate(double events, double exposure, double confidence = 0.95) {
  detail::require(events >= 0 && exposure > 0, "events ≥ 0 and exposure > 0");
  detail::require(confidence > 0 && confidence < 1, "a confidence in (0, 1)");
  const double alpha = 1 - confidence;
  Rate r;
  r.events = events, r.exposure = exposure, r.rate = events / exposure;
  r.lower = events > 0 ? tyche::probability::ChiSquared{2 * events}.quantile(alpha / 2) / (2 * exposure) : 0.0;
  r.upper = tyche::probability::ChiSquared{2 * events + 2}.quantile(1 - alpha / 2) / (2 * exposure);
  return r;
}

struct Standardised {
  double value = 0, se = 0, lower = 0, upper = 0;
  double observed = 0, expected = 0;  // indirect: D and E
};

// Direct standardisation of rates r_a = D_a / T_a with standard weights w_a.
inline Standardised direct(const std::vector<double>& events, const std::vector<double>& exposure, const std::vector<double>& weights, double confidence = 0.95) {
  detail::require(events.size() == exposure.size() && events.size() == weights.size() && !events.empty(), "events, exposure and a weight for each stratum");
  double W = 0, sum = 0, var = 0;
  for (std::size_t a = 0; a < events.size(); a++) {
    detail::require(exposure[a] > 0 && weights[a] >= 0, "exposure > 0 and weights ≥ 0");
    W += weights[a];
    sum += weights[a] * events[a] / exposure[a];
    var += weights[a] * weights[a] * events[a] / (exposure[a] * exposure[a]);
  }
  detail::require(W > 0, "weights with a positive sum");
  Standardised s;
  s.value = sum / W;
  s.se = std::sqrt(var) / W;
  const double z = tyche::probability::Normal{}.quantile(1 - (1 - confidence) / 2);
  s.lower = s.value - z * s.se, s.upper = s.value + z * s.se;
  return s;
}

// Indirect standardisation: observed events against those expected at the standard's rates.
inline Standardised indirect(double observed, const std::vector<double>& exposure, const std::vector<double>& standard_rates, double confidence = 0.95) {
  detail::require(exposure.size() == standard_rates.size() && !exposure.empty(), "an exposure and a standard rate for each stratum");
  double E = 0;
  for (std::size_t a = 0; a < exposure.size(); a++) E += exposure[a] * standard_rates[a];
  detail::require(E > 0, "events expected");
  const Rate r = rate(observed, E, confidence);
  Standardised s;
  s.observed = observed, s.expected = E, s.value = r.rate, s.lower = r.lower, s.upper = r.upper;
  s.se = std::sqrt(observed) / E;
  return s;
}

}  // namespace peira::cohorts
