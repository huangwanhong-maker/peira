// Exposure to treatments that reach beyond the units they are given to (peira's decision 0006).
//
// THE SETTING. Units i (sensors, places, plots, people) are observed in periods t. Treatments k begin in periods s_k and stay
// (absorbing). Each pair (i, k) has a distance d_ik ≥ 0 in whatever metric carries the spillover: metres along a street network,
// kilometres over the ground, steps in a graph. The caller computes the distances and gives them sparsely; a pair not given is
// beyond reach. The stable unit treatment value assumption (SUTVA; Rubin 1980) says a unit's outcome depends on its own treatment
// alone. Interference is its failure (Cox 1958; Hudgens and Halloran 2008, J. Am. Stat. Assoc. 103, 832–842).
//
// AN EXPOSURE MAPPING (Aronow and Samii 2017, Ann. Appl. Stat. 11, 1912–1947) reduces the whole assignment to a few conditions per
// unit. Here they are rings of distance with outer bounds b_0 < b_1 < … < b_{R−1}: unit i in period t is in ring r when the nearest
// treatment begun by t lies in [b_{r−1}, b_r) (b_{−1} = 0), and unexposed (ring R) when none lies nearer than b_{R−1}. Ring 0 is the
// direct exposure, rings 1 … R−1 the spillovers.
//
// The mapping is an assumption: the outcome depends on the assignment only through the ring and the period. In particular, nothing
// reaches beyond b_{R−1}, which is what makes the unexposed units comparisons. Within the mapping it cannot be tested. Finer rings
// show where estimated effects fall to zero (Butts 2021, §3).
//
// ONSETS. For each unit and ring r, the first period in which it is in ring r or nearer. onset(i, 0) begins its direct exposure, and
// onset(i, R − 1) is the first period it is reached at all. As treatments are absorbing, a unit's ring only comes nearer with time.
// Its entry ring is the ring in which it is first reached.

#pragma once

#include <peira/quasi/staggered.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace peira::interference {

using peira::quasi::never;

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::interference: ") + what);
}
}  // namespace detail

// How far a unit is from a treatment.
struct Contact {
  std::size_t unit = 0, treatment = 0;
  double distance = 0;
};

struct Exposure {
  std::vector<double> bounds;           // b_0 < … < b_{R−1}
  std::vector<std::int64_t> periods;    // the periods, ascending
  std::size_t units = 0;
  std::vector<int> ring;                // [i · T + t]: 0 … R − 1, or R where unexposed
  std::vector<std::int64_t> onset;      // [i · R + r]: the first period in ring r or nearer; never where none

  std::size_t rings() const { return bounds.size(); }
  int at(std::size_t i, std::size_t t) const { return ring[i * periods.size() + t]; }
  std::int64_t first(std::size_t i, std::size_t r) const { return onset[i * bounds.size() + r]; }
  std::int64_t reached(std::size_t i) const { return first(i, bounds.size() - 1); }
  // The ring in which unit i is first reached, or R if never.
  int entry(std::size_t i) const {
    const std::int64_t r0 = reached(i);
    if (r0 == never) return static_cast<int>(bounds.size());
    for (std::size_t r = 0; r < bounds.size(); r++)
      if (first(i, r) == r0) return static_cast<int>(r);
    return static_cast<int>(bounds.size());
  }
  // The index of a period's label.
  std::size_t index(std::int64_t period) const {
    const auto it = std::lower_bound(periods.begin(), periods.end(), period);
    detail::require(it != periods.end() && *it == period, "a period of the exposure");
    return static_cast<std::size_t>(it - periods.begin());
  }
};

// The ring of a distance under the bounds: the first r with d < b_r, or R.
inline int ring_of(double d, const std::vector<double>& bounds) {
  for (std::size_t r = 0; r < bounds.size(); r++)
    if (d < bounds[r]) return static_cast<int>(r);
  return static_cast<int>(bounds.size());
}

// The exposure of `units` units in the periods given (any order; kept ascending) to treatments beginning at `start[k]`, from the
// contacts given, under rings with the outer bounds given (ascending, positive).
inline Exposure expose(std::size_t units, std::vector<std::int64_t> periods, std::span<const std::int64_t> start, const std::vector<Contact>& contacts,
                       std::vector<double> bounds) {
  detail::require(!bounds.empty(), "one ring at least");
  for (std::size_t r = 0; r < bounds.size(); r++) detail::require(bounds[r] > 0 && (r == 0 || bounds[r] > bounds[r - 1]), "bounds positive and ascending");
  std::sort(periods.begin(), periods.end());
  periods.erase(std::unique(periods.begin(), periods.end()), periods.end());
  detail::require(!periods.empty(), "one period at least");
  Exposure e;
  e.bounds = std::move(bounds);
  e.periods = std::move(periods);
  e.units = units;
  const std::size_t T = e.periods.size(), R = e.bounds.size();
  // Per unit, its contacts as (start, ring).
  std::vector<std::vector<std::pair<std::int64_t, int>>> near(units);
  for (const auto& c : contacts) {
    detail::require(c.unit < units && c.treatment < start.size() && c.distance >= 0, "contacts of a unit, a treatment and a distance ≥ 0");
    const int r = ring_of(c.distance, e.bounds);
    if (r < static_cast<int>(R)) near[c.unit].push_back({start[c.treatment], r});
  }
  e.ring.assign(units * T, static_cast<int>(R));
  e.onset.assign(units * R, never);
  for (std::size_t i = 0; i < units; i++)
    for (std::size_t t = 0; t < T; t++) {
      int best = static_cast<int>(R);
      for (const auto& [s, r] : near[i])
        if (s <= e.periods[t]) best = std::min(best, r);
      e.ring[i * T + t] = best;
      for (int r = best; r < static_cast<int>(R); r++)
        if (e.onset[i * R + static_cast<std::size_t>(r)] == never) e.onset[i * R + static_cast<std::size_t>(r)] = e.periods[t];
    }
  return e;
}

// The ring of each observation of a long panel: `unit` gives each observation's unit as an index (0 … units − 1), `period` its label.
inline std::vector<int> rings_of(const Exposure& e, std::span<const std::size_t> unit, std::span<const std::int64_t> period) {
  detail::require(unit.size() == period.size(), "a unit and a period for every observation");
  std::vector<int> out(unit.size());
  for (std::size_t k = 0; k < unit.size(); k++) {
    detail::require(unit[k] < e.units, "units of the exposure");
    out[k] = e.at(unit[k], e.index(period[k]));
  }
  return out;
}

}  // namespace peira::interference
