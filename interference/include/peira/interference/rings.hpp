// Direct and spillover effects by rings of distance (peira's decision 0006), on the exposure of exposure.hpp.
//
// RINGS IN A TWO-WAY REGRESSION (Butts 2021, "Difference-in-differences estimation with spatial spillovers", arXiv:2105.03737;
// rings after Clarke 2017). With ring_it the ring of unit i in period t,
//   y_it = α_i + γ_t + Σ_{r < R} τ_r [ring_it = r] + ε_it,
// by the two-way within estimator with errors clustered by unit (quasi's two_way_fixed_effects, on tyche's panels). Three
// assumptions are needed:
// - parallel trends of the outcomes without treatment;
// - an exposure mapping that holds;
// - effects that are the same for every unit and period.
// Under them τ_0 is the direct effect and τ_r the spillover in ring r, each against the unexposed. With timing and effects that vary,
// staggered.hpp's caveat applies and ring_staggered below is the robust estimator.
//
// THE BIAS OF IGNORING SPILLOVERS. A regression on the direct indicator alone takes the units in the spillover rings as comparisons.
// In a 2 × 2 design with no noise its estimate is
//   τ_0 − Σ_{r ≥ 1} (n_r / n_c) τ_r,
// with n_r the units in ring r and n_c all the units it compares with (Butts 2021, §2). If the spillovers have the effect's sign,
// the effect is understated; if the opposite sign (trade drawn from neighbours), it is overstated. Checked exactly.
//
// RINGS WITH TIMING THAT VARIES (Callaway and Sant'Anna 2021, by ring; Butts 2021, §4). For ring r, each unit whose entry ring is r
// forms part of the group of its onset(i, r). It is compared with units not yet reached by any ring (or never reached). Its
// periods from the one in which it comes nearer than r are dropped, so that its effect is that of ring r alone. The group-time
// effects ATT_r(g, t) and their aggregations are quasi's staggered estimator with `reached`. Under parallel trends and the mapping,
// ATT_r is the average effect of entering ring r on those who enter it.
//
// SYNTHETIC CONTROL WITH CLEAN DONORS. A donor reached by a treatment carries its spillover into the synthetic unit and biases the
// gap. Cao and Dowd (2019, arXiv:1902.07343) model it. The remedy taken here leaves those donors out, as the mapping justifies:
// `clean_donors` returns the units never reached.

#pragma once

#include <peira/quasi/difference_in_differences.hpp>
#include <peira/quasi/staggered.hpp>

#include <calx/linalg/rectangular.hpp>
#include <calx/random/stream.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "exposure.hpp"

namespace peira::interference {

struct RingFit {
  std::vector<double> effect, se, t, p, lower, upper;  // per ring r < R: τ_r and its inference
  peira::quasi::WithinFit regression;
};

// The two-way regression on ring indicators: y and the unit and period labels of each observation, and its ring (rings_of).
inline RingFit ring_two_way(std::span<const double> y, std::span<const std::int64_t> unit, std::span<const std::int64_t> period, std::span<const int> ring,
                            std::size_t rings, double confidence = 0.95) {
  const std::size_t n = y.size();
  detail::require(unit.size() == n && period.size() == n && ring.size() == n, "a unit, a period and a ring for every observation");
  detail::require(rings >= 1, "one ring at least");
  std::vector<double> direct(n);
  calx::linalg::Rectangular others(n, rings - 1);
  for (std::size_t k = 0; k < n; k++) {
    detail::require(ring[k] >= 0 && ring[k] <= static_cast<int>(rings), "rings 0 … R");
    direct[k] = ring[k] == 0;
    for (std::size_t r = 1; r < rings; r++) others(k, r - 1) = ring[k] == static_cast<int>(r);
  }
  RingFit f;
  f.regression = peira::quasi::two_way_fixed_effects(y, direct, unit, period, others, confidence);
  const auto& L = f.regression.fit;
  for (std::size_t r = 0; r < rings; r++) {
    f.effect.push_back(L.beta[r]), f.se.push_back(L.se[r]), f.t.push_back(L.t[r]), f.p.push_back(L.p[r]);
    f.lower.push_back(L.ci_lower[r]), f.upper.push_back(L.ci_upper[r]);
  }
  return f;
}

// The staggered estimator for ring r: `unit` gives each observation's unit as an index of the exposure, `period` its label.
inline peira::quasi::StaggeredFit ring_staggered(std::span<const double> y, std::span<const std::size_t> unit, std::span<const std::int64_t> period, const Exposure& e,
                                                 std::size_t r, peira::quasi::Comparison comparison = peira::quasi::Comparison::not_yet_treated,
                                                 std::size_t replications = 0, calx::random::Stream* stream = nullptr) {
  const std::size_t n = y.size();
  detail::require(unit.size() == n && period.size() == n, "a unit and a period for every observation");
  detail::require(r < e.rings(), "a ring of the exposure");
  std::vector<double> ys;
  std::vector<std::int64_t> us, ps, first, reached;
  for (std::size_t k = 0; k < n; k++) {
    const std::size_t i = unit[k];
    detail::require(i < e.units, "units of the exposure");
    const bool target = e.entry(i) == static_cast<int>(r);
    // A unit of ring r leaves when it comes nearer.
    if (target && r > 0 && e.first(i, r - 1) != never && period[k] >= e.first(i, r - 1)) continue;
    ys.push_back(y[k]);
    us.push_back(static_cast<std::int64_t>(i));
    ps.push_back(period[k]);
    first.push_back(target ? e.first(i, r) : never);
    reached.push_back(e.reached(i));
  }
  return peira::quasi::staggered(ys, us, ps, first, comparison, replications, stream, reached);
}

// The units never reached by any treatment: the donors a synthetic control may use.
inline std::vector<std::size_t> clean_donors(const Exposure& e) {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < e.units; i++)
    if (e.reached(i) == never) out.push_back(i);
  return out;
}

}  // namespace peira::interference
