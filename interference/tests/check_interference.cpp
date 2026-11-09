// Checks of interference (exposure.hpp, rings.hpp, spillover_test.hpp).
//
// - The exposure by hand: rings, onsets and entry rings of four units under two treatments that begin at different times.
// - Rings in a two-way regression, on a panel without noise: the direct effect and the spillover recovered to rounding; the
//   regression that ignores spillovers off by exactly n_s/(n_s + n_u) τ_1; with noise, both within four standard errors.
// - Rings with timing that varies, without noise: ATT by ring recovered exactly against units not yet reached; the staggered
//   estimator that counts the spilled-over units as comparisons is biased.
// - Synthetic control: a donor touched by a spillover, the treated unit's twin before treatment, takes a weight w and biases the gap
//   by w times the spillover; with the clean donors alone the gap is the effect.
// - The test of no spillovers: exact by enumeration (P(p ≤ α) ≤ α in every class of the focal units' treatments, under the null);
//   and it rejects a spillover of 3 on a design of 24.

#include <peira/interference/interference.hpp>
#include <peira/quasi/difference_in_differences.hpp>
#include <peira/quasi/synthetic.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace pi = peira::interference;
namespace pq = peira::quasi;

namespace {

// A panel of units on a line at positions x_i, observed in periods 1 … T, with unit and period effects, treatments at sites s_k
// beginning at b_k, and constant effects by ring; noise of the sd given.
struct Panel {
  std::vector<double> y;
  std::vector<std::int64_t> unit_label, period;
  std::vector<std::size_t> unit;
  pi::Exposure e;
};

Panel simulate(const std::vector<double>& x, std::int64_t T, const std::vector<double>& site, const std::vector<std::int64_t>& begin, const std::vector<double>& bounds,
               const std::vector<double>& effect, double sd, calx::random::Stream& st) {
  std::vector<pi::Contact> contacts;
  for (std::size_t i = 0; i < x.size(); i++)
    for (std::size_t k = 0; k < site.size(); k++) contacts.push_back({i, k, std::fabs(x[i] - site[k])});
  std::vector<std::int64_t> periods;
  for (std::int64_t t = 1; t <= T; t++) periods.push_back(t);
  Panel p;
  p.e = pi::expose(x.size(), periods, begin, contacts, bounds);
  std::vector<double> alpha(x.size()), gamma(static_cast<std::size_t>(T));
  for (auto& a : alpha) a = st.uniform(0, 10);
  for (auto& g : gamma) g = st.uniform(-3, 3);
  for (std::size_t i = 0; i < x.size(); i++)
    for (std::int64_t t = 1; t <= T; t++) {
      const int r = p.e.at(i, static_cast<std::size_t>(t - 1));
      const double tau = r < static_cast<int>(effect.size()) ? effect[static_cast<std::size_t>(r)] : 0.0;
      p.y.push_back(alpha[i] + gamma[static_cast<std::size_t>(t - 1)] + tau + (sd > 0 ? st.normal(0, sd) : 0.0));
      p.unit_label.push_back(static_cast<std::int64_t>(i));
      p.unit.push_back(i);
      p.period.push_back(t);
    }
  return p;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  {
    // Four units, two treatments beginning in periods 2 and 4, rings of 100 and 500.
    const std::vector<pi::Contact> c{{0, 0, 50}, {1, 0, 300}, {1, 1, 80}, {2, 1, 450}};
    const auto e = pi::expose(4, {1, 2, 3, 4, 5}, std::vector<std::int64_t>{2, 4}, c, {100, 500});
    const std::vector<int> ring1{2, 1, 1, 0, 0};
    bool ok = true;
    for (std::size_t t = 0; t < 5; t++) ok = ok && e.at(1, t) == ring1[t] && e.at(3, t) == 2;
    s.require("exposure: a unit in the spillover ring from period 2 and directly exposed from 4; another never reached", ok);
    s.require("... its onsets 4 (ring 0) and 2 (ring 1), entering in ring 1; the first unit enters in ring 0 at 2; the third in ring 1 at 4",
              e.first(1, 0) == 4 && e.first(1, 1) == 2 && e.entry(1) == 1 && e.entry(0) == 0 && e.first(0, 0) == 2 && e.entry(2) == 1 && e.first(2, 0) == pi::never &&
                  e.entry(3) == 2);
  }
  {
    // Rings in a two-way regression: 120 units on a line, 4 sites treated from period 5 of 8, a direct effect of 5 within 100 m and a
    // spillover of 2 within 400 m.
    calx::random::Stream st(20261109, 0x696e7466u, 0, 0);  // "intf"
    std::vector<double> x;
    for (int i = 0; i < 120; i++) x.push_back(50.0 * i);
    const std::vector<double> site{400, 1500, 2600, 3700};
    const auto p = simulate(x, 8, site, {5, 5, 5, 5}, {100, 400}, {5, 2}, 0, st);
    const auto ring = pi::rings_of(p.e, p.unit, p.period);
    const auto f = pi::ring_two_way(p.y, p.unit_label, p.period, ring, 2);
    s.check("rings in a two-way regression, without noise: the direct effect 5 and the spillover 2", std::fabs(f.effect[0] - 5) + std::fabs(f.effect[1] - 2), 1e-9);
    // The regression on the direct indicator alone.
    std::vector<double> direct(p.y.size());
    for (std::size_t k = 0; k < ring.size(); k++) direct[k] = ring[k] == 0;
    const auto naive = pq::two_way_fixed_effects(p.y, direct, p.unit_label, p.period);
    std::size_t ns = 0, nu = 0;
    for (std::size_t i = 0; i < x.size(); i++) {
      const int r = p.e.at(i, 7);
      ns += r == 1, nu += r == 2;
    }
    const double expected = 5 - 2.0 * static_cast<double>(ns) / static_cast<double>(ns + nu);
    s.check("... ignoring spillovers takes 5 − 2 n_s/(n_s + n_u) (" + std::to_string(ns) + " units spilled over, " + std::to_string(nu) + " unexposed)",
            naive.fit.beta[0] - expected, 1e-9);
    const auto q = simulate(x, 8, site, {5, 5, 5, 5}, {100, 400}, {5, 2}, 1.0, st);
    const auto fq = pi::ring_two_way(q.y, q.unit_label, q.period, pi::rings_of(q.e, q.unit, q.period), 2);
    s.check("... with noise of sd 1, each within four standard errors (larger of the two |z|)",
            std::max(std::fabs(fq.effect[0] - 5) / fq.se[0], std::fabs(fq.effect[1] - 2) / fq.se[1]), 4);
  }
  {
    // Rings with timing that varies: sites treated from periods 4, 6 and 8 of 10; a direct effect of 5 and a spillover of 2.
    calx::random::Stream st(20261109, 0x73746167u, 0, 0);  // "stag"
    std::vector<double> x;
    for (int i = 0; i < 150; i++) x.push_back(20.0 * i);
    const std::vector<double> site{300, 1200, 2300};
    const auto p = simulate(x, 10, site, {4, 6, 8}, {100, 400}, {5, 2}, 0, st);
    const auto d = pi::ring_staggered(p.y, p.unit, p.period, p.e, 0);
    const auto sp = pi::ring_staggered(p.y, p.unit, p.period, p.e, 1);
    s.check("rings with timing that varies, without noise: ATT 5 in ring 0 and 2 in ring 1, against units not yet reached",
            std::fabs(d.simple - 5) + std::fabs(sp.simple - 2), 1e-9);
    double worst = 0;
    for (const auto& [ev, v] : d.by_exposure)
      if (ev < -1) worst = std::max(worst, std::fabs(v));
    s.check("... and the event study's placebo periods before entry are zero", worst, 1e-9);
    // The staggered estimator that takes the spilled-over units as never treated.
    std::vector<std::int64_t> first, label;
    for (std::size_t k = 0; k < p.unit.size(); k++) first.push_back(p.e.first(p.unit[k], 0)), label.push_back(static_cast<std::int64_t>(p.unit[k]));
    const auto plain = pq::staggered(p.y, label, p.period, first, pq::Comparison::not_yet_treated);
    s.require("... while counting the spilled-over units as comparisons misses the direct effect: " + std::to_string(plain.simple), std::fabs(plain.simple - 5) > 0.3);
  }
  {
    // Synthetic control: outcomes from two factors; the treated unit a convex combination of three clean donors, effect 4 from
    // period 9 of 14; a fourth donor its twin before period 9, touched by a spillover of 3 after.
    calx::random::Stream st(20261109, 0x73796e74u, 0, 0);  // "synt"
    const std::size_t T = 14, t0 = 8;
    std::vector<double> f1(T), f2(T);
    for (std::size_t t = 0; t < T; t++) f1[t] = st.uniform(5, 15), f2[t] = st.uniform(-4, 4);
    std::vector<std::vector<double>> donors;
    std::vector<std::pair<double, double>> load{{1.0, 0.2}, {0.4, 1.3}, {1.6, -0.8}, {0.9, 0.9}, {0.2, 0.1}};
    for (const auto& [a, b] : load) {
      std::vector<double> yj(T);
      for (std::size_t t = 0; t < T; t++) yj[t] = a * f1[t] + b * f2[t];
      donors.push_back(yj);
    }
    const double w[3] = {0.5, 0.3, 0.2};
    std::vector<double> treated(T), twin(T);
    for (std::size_t t = 0; t < T; t++) {
      const double base = w[0] * donors[0][t] + w[1] * donors[1][t] + w[2] * donors[2][t];
      treated[t] = base + (t >= t0 ? 4 : 0);
      twin[t] = base + (t >= t0 ? 3 : 0);
    }
    auto all = donors;
    all.push_back(twin);
    const auto biased = pq::synthetic_control(treated, all, t0, false);
    const auto clean = pq::synthetic_control(treated, donors, t0, false);
    double gb = 0, gc = 0;
    for (std::size_t t = t0; t < T; t++) gb += biased.gap[t] / static_cast<double>(T - t0), gc += clean.gap[t] / static_cast<double>(T - t0);
    const double wt = biased.weights.back();
    s.check("synthetic control: a donor spilled over, the treated unit's twin before, takes weight w and biases the gap to 4 − 3w (w = " + std::to_string(wt) + ")",
            gb - (4 - 3 * wt), 1e-4);
    s.require("... and w is not zero, so the bias is there", wt > 0.05);
    s.check("... and the clean donors alone show the effect, 4", gc - 4, 1e-4);
    // clean_donors on an exposure: the donors never reached.
    const auto e = pi::expose(6, {1, 2}, std::vector<std::int64_t>{2}, {{0, 0, 10}, {5, 0, 50}}, {100});
    s.require("... clean_donors leaves out the units a treatment reaches", pi::clean_donors(e) == std::vector<std::size_t>{1, 2, 3, 4});
  }
  {
    // The test of no spillovers, exact: 10 units on a line, neighbours adjacent, 3 treated; focal units the even ones; any outcomes.
    const std::size_t n = 10;
    std::vector<std::vector<std::size_t>> near(n);
    for (std::size_t i = 0; i < n; i++) {
      if (i > 0) near[i].push_back(i - 1);
      if (i + 1 < n) near[i].push_back(i + 1);
    }
    std::vector<int> focal(n);
    for (std::size_t i = 0; i < n; i++) focal[i] = i % 2 == 0;
    const std::vector<double> y{3.1, 0.4, 2.2, 5.0, 1.7, 2.9, 4.4, 0.8, 3.6, 2.5};
    std::map<std::vector<int>, std::vector<double>> classes;  // the focal units' treatments → the p-values of its assignments
    peira::design::each_complete(n, 3, [&](const peira::design::Assignment& z) {
      std::vector<int> key;
      for (std::size_t i = 0; i < n; i++)
        if (focal[i]) key.push_back(z[i]);
      classes[key].push_back(pi::no_spillover(y, z, focal, near).p_value);
    });
    double excess = 0;
    for (const auto& [key, ps] : classes)
      for (const double a : {0.05, 0.1, 0.2, 0.3, 0.5, 0.8}) {
        const double share = static_cast<double>(std::count_if(ps.begin(), ps.end(), [&](double p) { return p <= a + 1e-12; })) / static_cast<double>(ps.size());
        excess = std::max(excess, share - a);
      }
    s.check("the test of no spillovers is exact: P(p ≤ α) ≤ α in every class of the focal units' treatments (largest excess)", std::max(0.0, excess), 1e-12);
    // Power: 24 units on a line, 4 treated at 1, 7, 13 and 19, none of them focal (the focal units are the even ones, by a rule fixed
    // in advance); untreated neighbours of a treated unit raised by 3.
    const std::size_t m = 24;
    std::vector<std::vector<std::size_t>> nb(m);
    for (std::size_t i = 0; i < m; i++) {
      if (i > 0) nb[i].push_back(i - 1);
      if (i + 1 < m) nb[i].push_back(i + 1);
    }
    std::vector<int> z(m, 0), fm(m);
    for (const std::size_t i : {1, 7, 13, 19}) z[i] = 1;
    for (std::size_t i = 0; i < m; i++) fm[i] = i % 2 == 0;
    calx::random::Stream st(20261109, 0x706f7772u, 0, 0);  // "powr"
    std::vector<double> yy(m);
    for (std::size_t i = 0; i < m; i++) {
      bool ex = false;
      for (const std::size_t j : nb[i]) ex = ex || z[j] == 1;
      yy[i] = st.normal(0, 0.5) + (z[i] == 0 && ex ? 3.0 : 0.0);
    }
    const auto r = pi::no_spillover(yy, z, fm, nb);
    s.check("... and it finds a spillover of 3 on 24 units: p (exact, " + std::to_string(r.assignments) + " assignments)", r.p_value, 0.05);
  }
  return s.done();
}
