// An exact test of no spillovers in a randomised design (Athey, Eckles and Imbens 2018, "Exact p-values for network interference",
// J. Am. Stat. Assoc. 113, 230–240; peira's decision 0006).
//
// THE NULL. No unit's outcome depends on the treatments of the others. It is not sharp for every unit, since a treated unit's
// outcome under control is unseen. It is sharp for focal units F, chosen without looking at the outcomes, given their own
// treatment. Under the null, a focal unit's outcome is the same under every assignment that gives the focal units the treatments
// they received.
//
// THE TEST. The assignment is complete randomisation of n₁ treated among n units. Conditioning on the focal units' treatments
// leaves complete randomisation of the remaining treated among the non-focal units. Each such assignment is re-drawn, or
// enumerated when few. The statistic on F is recomputed with the focal outcomes held fixed, and the p-value is the share of
// assignments whose statistic is at least as extreme as the one observed. It is exact: under the null P(p ≤ α) ≤ α for every α, in
// every conditioning class (checked by enumeration).
//
// THE STATISTIC. The difference in mean outcome between the untreated focal units that are exposed (near a treated unit, by the
// caller's relation: within a distance, an edge) and those that are not. It is AEI's simplest contrast. Where one group is empty,
// the statistic is 0. Any function of the focal outcomes and the assignment keeps the test valid; this one has power against
// spillovers on the untreated.

#pragma once

#include <peira/design/design.hpp>

#include <calx/random/stream.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "exposure.hpp"

namespace peira::interference {

struct SpilloverTest {
  double observed = 0;          // the statistic on the assignment observed
  double p_value = 1;           // two-sided
  std::size_t assignments = 0;  // enumerated, or drawn
  bool exact = false;
};

namespace detail {
// The difference in mean outcome between untreated focal units near a treated unit and those not near one.
inline double contrast(std::span<const double> y, std::span<const int> z, std::span<const int> focal, const std::vector<std::vector<std::size_t>>& near) {
  double s1 = 0, s0 = 0;
  std::size_t n1 = 0, n0 = 0;
  for (std::size_t i = 0; i < y.size(); i++) {
    if (!focal[i] || z[i]) continue;
    bool exposed = false;
    for (const std::size_t j : near[i]) exposed = exposed || z[j] == 1;
    if (exposed) s1 += y[i], n1++;
    else s0 += y[i], n0++;
  }
  return (n1 > 0 && n0 > 0) ? s1 / static_cast<double>(n1) - s0 / static_cast<double>(n0) : 0.0;
}
}  // namespace detail

// The test, for outcomes y, the assignment z observed (complete randomisation), the focal units (1 for focal) and, for each unit,
// the units near it. Exact by enumeration when the conditional design has at most `max_exact` assignments, else `draws` draws
// from the stream given.
inline SpilloverTest no_spillover(std::span<const double> y, std::span<const int> z, std::span<const int> focal, const std::vector<std::vector<std::size_t>>& near,
                                  std::size_t max_exact = 200000, std::size_t draws = 20000, calx::random::Stream* stream = nullptr) {
  const std::size_t n = y.size();
  detail::require(z.size() == n && focal.size() == n && near.size() == n, "an assignment, a focal mark and neighbours for every unit");
  for (const auto& v : near)
    for (const std::size_t j : v) detail::require(j < n, "neighbours among the units");
  // The non-focal units, among which the treated not focal are re-drawn.
  std::vector<std::size_t> free;
  std::size_t k = 0;
  for (std::size_t i = 0; i < n; i++)
    if (!focal[i]) free.push_back(i), k += z[i] == 1;
  SpilloverTest r;
  r.observed = detail::contrast(y, z, focal, near);
  const double tol = 1e-12 * (1 + std::fabs(r.observed));
  std::vector<int> w(z.begin(), z.end());
  std::size_t total = 0, hits = 0;
  const auto score = [&](const peira::design::Assignment& sub) {
    for (std::size_t m = 0; m < free.size(); m++) w[free[m]] = sub[m];
    total++;
    hits += std::fabs(detail::contrast(y, w, focal, near)) >= std::fabs(r.observed) - tol;
  };
  if (peira::design::binomial_coefficient(free.size(), k) <= static_cast<double>(max_exact)) {
    peira::design::each_complete(free.size(), k, score);
    r.exact = true;
  } else {
    detail::require(stream != nullptr, "a stream when the design is too large to enumerate");
    for (std::size_t b = 0; b < draws; b++) score(peira::design::complete(free.size(), k, *stream));
  }
  r.assignments = total;
  r.p_value = static_cast<double>(hits) / static_cast<double>(total);
  return r;
}

}  // namespace peira::interference
