// Checks of the quasi-experiments (difference_in_differences.hpp, staggered.hpp, instrumental.hpp's Wald estimator,
// discontinuity.hpp, synthetic.hpp).
//
// - Differences in differences (moved from nomos): the 2 × 2 estimate is the difference of the four means and recovers the effect;
//   two-way fixed effects on the balanced 2 × 2 panel give the same estimate.
// - Staggered adoption, noise-free, with effects that differ by group and grow with exposure: every ATT(g, t) exactly, the
//   placebos before treatment zero, the aggregation by exposure the size-weighted mean of the true effects, with never-treated or
//   not-yet-treated comparisons; the static two-way regression departs from the average effect, while with one constant effect it
//   agrees exactly; with homogeneous dynamics the two-way event study recovers each θ(e); with noise, the bootstrap's errors cover.
// - The Wald estimator is the ratio of the differences of means exactly, and recovers the compliers' effect among always-takers
//   and never-takers with other effects.
// - Regression discontinuity: lines on each side give the jump exactly at any bandwidth and kernel; curvature gives a bias that falls
//   fourfold when the bandwidth halves; with noise the effect lies within four standard errors; the fuzzy estimate is the ratio of
//   the two jumps and recovers the effect exactly when the outcome is linear in the treatment; cross-validation picks a narrower
//   bandwidth for a curved mean than for a straight one.
// - Synthetic control: a treated unit that is a convex combination of donors before treatment gets those weights back, and its gaps
//   are the effect; the placebos rank it first, p = 1/(J + 1); on any target the weights satisfy the optimality conditions on the
//   simplex; the projection onto the simplex.

#include <peira/quasi/difference_in_differences.hpp>
#include <peira/quasi/discontinuity.hpp>
#include <peira/quasi/instrumental.hpp>
#include <peira/quasi/staggered.hpp>
#include <peira/quasi/synthetic.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace pq = peira::quasi;
using calx::verification::exponential;

namespace {
calx::random::Stream z(20261108, 0x71756173u, 0, 0);  // "quas"

struct Panel {
  std::vector<double> y;
  std::vector<std::int64_t> unit, period, first;
};

// Units in groups first treated at 4 and 7, and never; periods 1..10; y = α_i + λ_t + τ_g(t − g) for t ≥ g, plus noise.
Panel staggered_panel(const std::function<double(std::int64_t, std::int64_t)>& tau, double noise, std::size_t per_group = 60) {
  Panel p;
  const std::int64_t groups[3] = {4, 7, pq::never};
  std::vector<double> lambda(11);
  for (auto& l : lambda) l = z.normal();
  std::int64_t id = 0;
  for (const std::int64_t g : groups)
    for (std::size_t k = 0; k < per_group; k++, id++) {
      const double a = 2 * z.normal();
      for (std::int64_t t = 1; t <= 10; t++) {
        double y = a + lambda[static_cast<std::size_t>(t)];
        if (g != pq::never && t >= g) y += tau(g, t - g);
        y += noise * z.normal();
        p.y.push_back(y), p.unit.push_back(id), p.period.push_back(t), p.first.push_back(g);
      }
    }
  return p;
}
}  // namespace

int main() {
  calx::verification::Suite s;
  // ---- differences in differences, 2 × 2 (moved from nomos's check_panel) ------------------------------------------
  {
    const std::size_t N = 120;
    std::vector<double> y, D;
    std::vector<int> treated, post;
    std::vector<std::int64_t> unit, period;
    for (std::size_t i = 0; i < N; i++) {
      const int g = i < N / 2 ? 1 : 0;
      const double a = z.normal();
      for (int t = 0; t < 2; t++) {
        y.push_back(a + 0.3 * g + 0.5 * t + 1.2 * g * t + 0.4 * z.normal());
        treated.push_back(g), post.push_back(t), D.push_back(g * t);
        unit.push_back(static_cast<std::int64_t>(i)), period.push_back(t);
      }
    }
    const auto did = pq::difference_in_differences(y, treated, post, unit);
    const double four = (did.mean[1][1] - did.mean[1][0]) - (did.mean[0][1] - did.mean[0][0]);
    s.check("differences in differences, 2 × 2: the estimate is the difference of the four means", did.estimate - four, 1e-12);
    s.require("... and the effect 1.2 lies within four standard errors", std::fabs(did.estimate - 1.2) < 4 * did.se, exponential(did.estimate) + " ± " + exponential(did.se));
    const auto twfe = pq::two_way_fixed_effects(y, D, unit, period);
    s.check("... two-way fixed effects on the balanced 2 × 2 panel give the same estimate", twfe.fit.beta[0] - did.estimate, 1e-10);
  }
  // ---- staggered adoption ------------------------------------------------------------------------------------------
  {
    const auto tau = [](std::int64_t g, std::int64_t e) { return g == 4 ? 1.0 + 0.5 * static_cast<double>(e) : 2.0 + 0.2 * static_cast<double>(e); };
    const auto P = staggered_panel(tau, 0.0);
    for (const auto comparison : {pq::Comparison::never_treated, pq::Comparison::not_yet_treated}) {
      const auto f = pq::staggered(P.y, P.unit, P.period, P.first, comparison);
      double post = 0, pre = 0;
      std::size_t npost = 0, npre = 0;
      for (const auto& e : f.effects) {
        if (e.period >= e.group) post = std::max(post, std::fabs(e.att - tau(e.group, e.period - e.group))), npost++;
        else pre = std::max(pre, std::fabs(e.att)), npre++;
      }
      const std::string which = comparison == pq::Comparison::never_treated ? "never treated" : "not yet treated";
      s.check("staggered, noise-free, against the " + which + ": every ATT(g, t) after treatment is the true effect (" + std::to_string(npost) + ")", post, 1e-10);
      s.check("... and every placebo before treatment is zero (" + std::to_string(npre) + ")", pre, 1e-10);
    }
    const auto f = pq::staggered(P.y, P.unit, P.period, P.first);
    double worst = 0;
    for (const auto& [e, v] : f.by_exposure) {
      double num = 0, den = 0;
      for (const std::int64_t g : {4, 7})
        if (g + e >= 1 && g + e <= 10 && g + e != g - 1) num += 60 * (e >= 0 ? tau(g, e) : 0.0), den += 60;
      worst = std::max(worst, std::fabs(v - num / den));
    }
    s.check("... the aggregation by exposure is the size-weighted mean of the true effects", worst, 1e-10);
    double truth = 0, n = 0;
    for (const std::int64_t g : {4, 7})
      for (std::int64_t t = g; t <= 10; t++) truth += tau(g, t - g), n += 1;
    s.check("... and the simple aggregation the mean of the true effects after treatment (groups of equal size)", f.simple - truth / n, 1e-10);
    const auto tw = pq::staggered_two_way(P.y, P.unit, P.period, P.first);
    s.require("... the static two-way regression departs from it by more than a tenth (dynamic effects, staggered timing)",
              std::fabs(tw.fit.beta[0] - f.simple) > 0.1 * std::fabs(f.simple), exponential(tw.fit.beta[0]) + " against " + exponential(f.simple));
    const auto C = staggered_panel([](std::int64_t, std::int64_t) { return 1.5; }, 0.0);
    const auto fc = pq::staggered(C.y, C.unit, C.period, C.first);
    const auto twc = pq::staggered_two_way(C.y, C.unit, C.period, C.first);
    s.check("with one constant effect, the two-way regression and the simple aggregation are both 1.5", std::fabs(twc.fit.beta[0] - 1.5) + std::fabs(fc.simple - 1.5), 1e-9);
    const auto H = staggered_panel([](std::int64_t, std::int64_t e) { return 0.5 + 0.25 * static_cast<double>(e); }, 0.0);
    const auto es = pq::event_study_two_way(H.y, H.unit, H.period, H.first, -6, 6);
    double ew = 0;
    for (const auto& [e, b] : es.coefficient) ew = std::max(ew, std::fabs(b - (e >= 0 ? 0.5 + 0.25 * static_cast<double>(e) : 0.0)));
    s.check("with dynamics the same in every group, the two-way event study recovers each θ(e) (leads zero)", ew, 1e-8);
    const auto N = staggered_panel(tau, 1.0, 80);
    calx::random::Stream boot(7, 0x626f6f74u, 0, 0);
    const auto fn = pq::staggered(N.y, N.unit, N.period, N.first, pq::Comparison::never_treated, 200, &boot);
    double tn = 0, nn = 0;
    for (const std::int64_t g : {4, 7})
      for (std::int64_t t = g; t <= 10; t++) tn += tau(g, t - g), nn += 1;
    s.require("with noise, the simple aggregation lies within four bootstrap standard errors of the truth",
              fn.simple_se > 0 && std::fabs(fn.simple - tn / nn) < 4 * fn.simple_se, exponential(fn.simple) + " ± " + exponential(fn.simple_se));
  }
  // ---- the Wald estimator -----------------------------------------------------------------------------------------
  {
    const std::size_t n = 20000;
    std::vector<double> y(n), d(n);
    std::vector<int> zi(n);
    for (std::size_t i = 0; i < n; i++) {
      const double u = z.uniform();
      const int type = u < 0.2 ? 0 : u < 0.5 ? 1 : 2;  // always-taker, never-taker, complier
      zi[i] = z.uniform() < 0.5;
      d[i] = type == 0 ? 1 : type == 1 ? 0 : zi[i];
      const double effect = type == 2 ? 2.0 : 5.0;
      y[i] = 1 + 0.5 * type + effect * d[i] + z.normal();
    }
    const auto w = pq::wald(y, d, zi);
    s.check("Wald: 2SLS with a binary instrument is the ratio of the differences of means", w.estimate - w.reduced_form / w.compliers, 1e-10);
    s.require("... and recovers the compliers' effect 2 (always-takers' is 5) within four standard errors", std::fabs(w.estimate - 2.0) < 4 * w.se,
              exponential(w.estimate) + " ± " + exponential(w.se) + ", compliers " + exponential(w.compliers));
  }
  // ---- regression discontinuity -----------------------------------------------------------------------------------
  {
    std::vector<double> x, y;
    for (int i = 0; i <= 2000; i++) {
      const double v = -1 + i * 0.001;
      x.push_back(v);
      const double D = v >= 0 ? 1.0 : 0.0;
      y.push_back(1 + 0.5 * v + 2 * D + 0.3 * D * v);
    }
    double worst = 0;
    for (const auto k : {pq::Kernel::triangular, pq::Kernel::epanechnikov, pq::Kernel::uniform})
      for (const double h : {0.2, 0.5}) worst = std::max(worst, std::fabs(pq::sharp_discontinuity(x, y, 0, h, k).estimate - 2));
    s.check("discontinuity: lines on each side give the jump 2 exactly, at every bandwidth and kernel", worst, 1e-10);
    std::vector<double> q(x.size());
    for (std::size_t i = 0; i < x.size(); i++) q[i] = 1 + (x[i] >= 0 ? 1.0 + x[i] * x[i] : -x[i] * x[i]);
    const double b1 = pq::sharp_discontinuity(x, q, 0, 0.4).estimate - 1, b2 = pq::sharp_discontinuity(x, q, 0, 0.2).estimate - 1;
    s.check("... curvature gives a bias that falls fourfold when the bandwidth halves (ratio − 4)", b1 / b2 - 4, 0.02);
    std::vector<double> xn, yn;
    for (int i = 0; i < 4000; i++) {
      const double v = z.uniform(-1, 1);
      xn.push_back(v), yn.push_back(1 + v + 0.8 * (v >= 0) + 0.5 * z.normal());
    }
    const auto r = pq::sharp_discontinuity(xn, yn, 0, 0.3);
    s.require("... with noise, the effect 0.8 lies within four standard errors", std::fabs(r.estimate - 0.8) < 4 * r.se, exponential(r.estimate) + " ± " + exponential(r.se));
    std::vector<double> xf, yf, df;
    for (int i = 0; i < 3000; i++) {
      const double v = z.uniform(-1, 1), D = v >= 0 ? 1.0 : 0.0;
      const double take = z.uniform() < 0.3 + 0.4 * D ? 1.0 : 0.0;
      xf.push_back(v), df.push_back(take), yf.push_back(1 + 3 * take + 0.5 * v + 0.2 * D * v);
    }
    const auto f = pq::fuzzy_discontinuity(xf, yf, df, 0, 0.5);
    s.check("fuzzy: the estimate is the ratio of the two local jumps", f.estimate - f.jump_outcome / f.jump_treatment, 1e-10);
    s.check("... and recovers the effect 3 exactly when the outcome is linear in the treatment", f.estimate - 3, 1e-9);
    std::vector<double> xc, yc, yl;
    for (int i = 0; i < 1500; i++) {
      const double v = z.uniform(-1, 1), e = 0.2 * z.normal();
      xc.push_back(v), yc.push_back(std::sin(4 * v) + (v >= 0) + e), yl.push_back(0.5 * v + (v >= 0) + e);
    }
    const std::vector<double> grid{0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8};
    const auto bc = pq::cross_validated_bandwidth(xc, yc, 0, grid), bl = pq::cross_validated_bandwidth(xc, yl, 0, grid);
    s.require("cross-validation picks a narrower bandwidth for a curved mean than for a straight one", bc.chosen < bl.chosen,
              std::to_string(bc.chosen) + " against " + std::to_string(bl.chosen));
  }
  // ---- synthetic control ------------------------------------------------------------------------------------------
  {
    const std::size_t J = 5, T = 20, T0 = 12;
    std::vector<std::vector<double>> donors(J, std::vector<double>(T));
    for (auto& d : donors) {
      double v = z.normal();
      for (auto& x : d) x = (v += z.normal());
    }
    std::vector<double> treated(T);
    for (std::size_t t = 0; t < T; t++) treated[t] = 0.5 * donors[0][t] + 0.3 * donors[2][t] + 0.2 * donors[3][t] + (t >= T0 ? 2.0 : 0.0);
    const auto sc = pq::synthetic_control(treated, donors, T0);
    const std::vector<double> truth{0.5, 0, 0.3, 0.2, 0};
    double dw = 0, dg = 0;
    for (std::size_t j = 0; j < J; j++) dw = std::max(dw, std::fabs(sc.weights[j] - truth[j]));
    for (std::size_t t = T0; t < T; t++) dg = std::max(dg, std::fabs(sc.gap[t] - 2.0));
    s.check("synthetic control: a convex combination of donors gets its weights back", dw, 1e-6);
    s.check("... its gaps after treatment are the effect 2", dg, 1e-5);
    s.require("... and the placebos rank it first: p = 1/(J + 1)", std::fabs(sc.p_value - 1.0 / 6) < 1e-12, std::to_string(sc.p_value));
    std::vector<double> other(T);
    for (auto& v : other) v = 3 * z.normal();
    const auto so = pq::synthetic_control(other, donors, T0, false);
    std::vector<double> g(J, 0.0);
    for (std::size_t j = 0; j < J; j++)
      for (std::size_t t = 0; t < T0; t++) g[j] += -2 * donors[j][t] * so.gap[t];
    const double gmin = *std::min_element(g.begin(), g.end());
    double kkt = 0;
    for (std::size_t j = 0; j < J; j++) kkt = std::max(kkt, so.weights[j] > 1e-8 ? std::fabs(g[j] - gmin) : std::max(0.0, gmin - g[j]));
    s.check("... on any target, the weights satisfy the simplex's optimality conditions (the gradient least on the support)", kkt / (1 + std::fabs(gmin)), 1e-6);
    const auto p = pq::project_simplex({0.3, -0.2, 1.5, 0.1});
    const auto same = pq::project_simplex({0.2, 0.3, 0.5});
    s.check("the projection onto the simplex sums to one and fixes the simplex's points",
            std::fabs(std::accumulate(p.begin(), p.end(), 0.0) - 1) + std::fabs(same[0] - 0.2) + std::fabs(same[1] - 0.3) + std::fabs(same[2] - 0.5), 1e-12);
  }
  return s.done();
}
