// Checks of peira/detection: closed populations by capture and recapture (capture_recapture.hpp).
//
// - Lincoln–Petersen and Chapman by hand (n1 = 50, n2 = 60, m2 = 15): 200 with variance 1400, and 3095/16 with Seber's variance
//   4899825/4352; Chao's interval satisfies (L − M)(U − M) = f̂0², and holds N̂ within it and above M.
// - Chapman's estimator is exactly unbiased under the hypergeometric law when n1 + n2 ≥ N, by enumerating m2; when n1 + n2 < N
//   its bias is exactly −(n2 + 1) C(N − n1, n2 + 1)/C(N, n2) (docs/theory.md, Theorem 10.5).
// - Schnabel by hand on four occasions, and equal to Lincoln–Petersen on two; the summary of capture histories by hand.
// - Models M0 and Mt: the profile likelihood equals the likelihood maximised over p on a grid; N̂ and the profile interval
//   equal the exhaustive search over every integer N from M to 50M; with no recapture N̂ = ∞.
// - Recovery from 300 simulated populations of 300 (five occasions): N̂'s bias below a quarter of its spread and the profile
//   intervals' coverage within 3.5 Monte Carlo standard errors of 95%, for M0 on constant p and Mt on varying p; AIC chooses Mt
//   when p varies and M0 mostly when it does not.

#include <peira/detection/detection.hpp>

#include <tyche/probability/distributions.hpp>

#include <calx/random/stream.hpp>
#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace pd = peira::detection;

namespace {

// E[Chapman's N̂] under m2 ~ Hypergeometric(N, n1, n2), by enumeration.
double chapman_expectation(int N, int n1, int n2) {
  double e = 0;
  for (int m = std::max(0, n1 + n2 - N); m <= std::min(n1, n2); m++) {
    const double lp = calx::special_functions::log_choose(n1, m) + calx::special_functions::log_choose(N - n1, n2 - m) - calx::special_functions::log_choose(N, n2);
    e += std::exp(lp) * ((n1 + 1.0) * (n2 + 1.0) / (m + 1.0) - 1);
  }
  return e;
}

struct Exhaustive {
  double argmax = 0, lower = 0, upper = 0;
};

// Every integer N in [M, top]: the profile's maximiser and its interval.
template <class F>
Exhaustive exhaustive(F&& profile, double M, double top, double confidence) {
  const double cut = tyche::probability::ChiSquared{1}.upper_quantile(1 - confidence) / 2;
  std::vector<double> v;
  double best = -std::numeric_limits<double>::infinity();
  Exhaustive e;
  for (double N = M; N <= top; N++) {
    v.push_back(profile(N));
    if (v.back() > best) best = v.back(), e.argmax = N;
  }
  e.lower = top, e.upper = M;
  for (std::size_t k = 0; k < v.size(); k++)
    if (best - v[k] <= cut) e.lower = std::min(e.lower, M + k), e.upper = std::max(e.upper, M + k);
  return e;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Two samples.
  {
    const auto lp = pd::lincoln_petersen(50, 60, 15), ch = pd::chapman(50, 60, 15);
    s.check("lincoln-petersen: N̂ = n1 n2/m2 = 200", lp.abundance - 200, 1e-12);
    s.check("lincoln-petersen: variance n1 n2 (n1 − m2)(n2 − m2)/m2³ = 1400", lp.se * lp.se - 1400, 1e-9);
    s.check("chapman: N̂ = (n1 + 1)(n2 + 1)/(m2 + 1) − 1 = 3095/16", ch.abundance - 3095.0 / 16, 1e-12);
    s.check("chapman: Seber's variance = 4899825/4352", ch.se * ch.se - 4899825.0 / 4352, 1e-9);
    const double f0 = ch.abundance - ch.seen;
    s.check("chao's interval: (L − M)(U − M) = f̂0²", (ch.lower - ch.seen) * (ch.upper - ch.seen) / (f0 * f0) - 1, 1e-12);
    s.require("chao's interval: M < L < N̂ < U", ch.seen < ch.lower && ch.lower < ch.abundance && ch.abundance < ch.upper,
              "[" + std::to_string(ch.lower) + ", " + std::to_string(ch.upper) + "]");
    double worst = 0;
    for (const auto& c : std::vector<std::array<int, 3>>{{30, 12, 18}, {30, 15, 15}, {30, 20, 15}, {50, 30, 25}, {80, 40, 41}})
      worst = std::max(worst, std::fabs(chapman_expectation(c[0], c[1], c[2]) - c[0]));
    s.check("chapman: exactly unbiased when n1 + n2 ≥ N (five cases, by enumeration)", worst, 1e-10);
    double worst_bias = 0;
    for (const auto& c : std::vector<std::array<int, 3>>{{30, 10, 15}, {40, 8, 12}, {60, 20, 25}}) {
      const int N = c[0], n1 = c[1], n2 = c[2];
      const double bias = -(n2 + 1) * std::exp(calx::special_functions::log_choose(N - n1, n2 + 1) - calx::special_functions::log_choose(N, n2));
      worst_bias = std::max(worst_bias, std::fabs(chapman_expectation(N, n1, n2) - N - bias));
    }
    s.check("chapman: when n1 + n2 < N the bias is −(n2 + 1) C(N − n1, n2 + 1)/C(N, n2) (three cases, from 10⁻⁵ to 0.8)", worst_bias, 1e-10);
  }
  // ---- Schnabel, and summaries of histories.
  {
    const std::vector<double> C{10, 15, 12, 20}, R{0, 3, 4, 9};
    const auto sc = pd::schnabel(C, R);
    s.check("schnabel: Σ C_t M_t / Σ R_t = 1014/16 by hand", sc.abundance - 1014.0 / 16, 1e-12);
    s.check("schnabel: Chapman's form Σ C_t M_t / (Σ R_t + 1) = 1014/17", sc.abundance_chapman - 1014.0 / 17, 1e-12);
    s.require("schnabel: the Poisson interval holds N̂", sc.lower < sc.abundance && sc.abundance < sc.upper,
              "[" + std::to_string(sc.lower) + ", " + std::to_string(sc.upper) + "]");
    const auto two = pd::schnabel(std::vector<double>{50, 60}, std::vector<double>{0, 15});
    s.check("schnabel: on two occasions it is Lincoln–Petersen", two.abundance - 200, 1e-12);
    pd::CaptureHistories h(4, 3);
    const int rows[4][3] = {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {0, 0, 1}};
    for (std::size_t i = 0; i < 4; i++)
      for (std::size_t j = 0; j < 3; j++) h(i, j) = rows[i][j];
    const auto sum = pd::summarise(h);
    const bool ok = sum.caught == std::vector<double>{2, 2, 3} && sum.recaptured == std::vector<double>{0, 1, 2} &&
                    sum.marked_before == std::vector<double>{0, 2, 3} && sum.distinct == 4 && sum.captures == 7;
    s.require("summarise: caught (2, 2, 3), recaptured (0, 1, 2), marked before (0, 2, 3), 4 distinct, 7 captures", ok);
  }
  // ---- M0 and Mt against exhaustive searches.
  {
    calx::random::Stream st(1978, 1, 0, 0);
    const std::vector<double> pt{0.15, 0.3, 0.2, 0.35, 0.25};
    const auto h = pd::simulate_captures(st, 120, pt);
    const auto sm = pd::summarise(h);
    const double M = sm.distinct, t = 5;
    const auto m0 = pd::closed_m0(h), mt = pd::closed_mt(h);
    auto prof0 = [&](double N) { std::vector<double> p(5, sm.captures / (t * N)); return pd::closed_log_likelihood(N, p, sm.caught, M); };
    auto proft = [&](double N) {
      std::vector<double> p(5);
      for (std::size_t j = 0; j < 5; j++) p[j] = sm.caught[j] / N;
      return pd::closed_log_likelihood(N, p, sm.caught, M);
    };
    const auto e0 = exhaustive(prof0, M, 50 * M, 0.95), et = exhaustive(proft, M, 50 * M, 0.95);
    s.check("M0: N̂ equals the exhaustive search over N in [M, 50M]", m0.abundance - e0.argmax, 0);
    s.check("M0: the profile interval equals the exhaustive one", std::fabs(m0.lower - e0.lower) + std::fabs(m0.upper - e0.upper), 0);
    s.check("Mt: N̂ equals the exhaustive search over N in [M, 50M]", mt.abundance - et.argmax, 0);
    s.check("Mt: the profile interval equals the exhaustive one", std::fabs(mt.lower - et.lower) + std::fabs(mt.upper - et.upper), 0);
    s.check("M0: the log-likelihood reported is the profile's at N̂", m0.log_likelihood - prof0(m0.abundance), 1e-12);
    // p̂ = n/(tN) maximises the likelihood over p at N̂: a grid of 10⁵ values of p.
    double grid0 = -std::numeric_limits<double>::infinity(), gridt = 0;
    for (int k = 1; k < 100000; k++) {
      const double p = k / 100000.0;
      grid0 = std::max(grid0, pd::closed_log_likelihood(m0.abundance, std::vector<double>(5, p), sm.caught, M));
    }
    gridt = calx::special_functions::log_gamma(mt.abundance + 1) - calx::special_functions::log_gamma(mt.abundance - M + 1);
    for (std::size_t j = 0; j < 5; j++) {
      double b = -std::numeric_limits<double>::infinity();
      for (int k = 1; k < 100000; k++) {
        const double p = k / 100000.0;
        b = std::max(b, sm.caught[j] * std::log(p) + (mt.abundance - sm.caught[j]) * std::log1p(-p));
      }
      gridt += b;
    }
    s.check("M0: the profile at N̂ equals the likelihood maximised over p on a grid (≥, to the grid's resolution)", m0.log_likelihood - grid0, 1e-6);
    s.check("Mt: the profile at N̂ equals the likelihood maximised over each p_j on a grid", mt.log_likelihood - gridt, 1e-6);
    s.require("M0, Mt: the grid never exceeds the profile", grid0 <= m0.log_likelihood + 1e-12 && gridt <= mt.log_likelihood + 1e-12);
    s.check("Mt: the probability of being seen is 1 − Π(1 − n_j/N̂)", mt.probability_seen - (1 - (1 - sm.caught[0] / mt.abundance) * (1 - sm.caught[1] / mt.abundance) *
                                                                                              (1 - sm.caught[2] / mt.abundance) * (1 - sm.caught[3] / mt.abundance) *
                                                                                              (1 - sm.caught[4] / mt.abundance)),
            1e-14);
    // No recapture: nobody caught twice.
    const auto none = pd::closed_m0(10, 10, 3, 0.95, 5000);
    s.require("M0: with no recapture N̂ and the upper limit are ∞", std::isinf(none.abundance) && std::isinf(none.upper));
  }
  // ---- Recovery from simulated populations.
  {
    const std::size_t R = 300;
    const double N = 300;
    const std::vector<double> pc(5, 0.2), pv{0.1, 0.35, 0.15, 0.4, 0.2};
    double s0 = 0, s02 = 0, c0 = 0, st_ = 0, st2 = 0, ct = 0, mt_chosen = 0, m0_chosen = 0;
    for (std::size_t rep = 0; rep < R; rep++) {
      calx::random::Stream a(1978, 2, static_cast<std::uint32_t>(rep), 0), b(1978, 3, static_cast<std::uint32_t>(rep), 0);
      const auto hc = pd::simulate_captures(a, static_cast<std::size_t>(N), pc), hv = pd::simulate_captures(b, static_cast<std::size_t>(N), pv);
      const auto f0 = pd::closed_m0(hc), ft = pd::closed_mt(hv);
      s0 += f0.abundance, s02 += f0.abundance * f0.abundance, c0 += f0.lower <= N && N <= f0.upper;
      st_ += ft.abundance, st2 += ft.abundance * ft.abundance, ct += ft.lower <= N && N <= ft.upper;
      mt_chosen += ft.aic < pd::closed_m0(hv).aic;
      m0_chosen += f0.aic < pd::closed_mt(hc).aic;
    }
    const double n = static_cast<double>(R), lim = 3.5 * std::sqrt(0.95 * 0.05 / n);
    const double m0 = s0 / n, sd0 = std::sqrt(s02 / n - m0 * m0), mt = st_ / n, sdt = std::sqrt(st2 / n - mt * mt);
    s.check("M0 recovery: N̂'s bias / its spread (Monte Carlo error 0.058)", (m0 - N) / sd0, 0.25);
    s.check("M0 recovery: the profile interval's coverage − 0.95 (limit 3.5 Monte Carlo s.e.)", c0 / n - 0.95, lim);
    s.check("Mt recovery: N̂'s bias / its spread (Monte Carlo error 0.058)", (mt - N) / sdt, 0.25);
    s.check("Mt recovery: the profile interval's coverage − 0.95 (limit 3.5 Monte Carlo s.e.)", ct / n - 0.95, lim);
    s.require("AIC: Mt chosen in at least 95% of populations whose p varies", mt_chosen / n >= 0.95, std::to_string(mt_chosen / n));
    s.require("AIC: M0 chosen in at least 80% of populations whose p does not", m0_chosen / n >= 0.8, std::to_string(m0_chosen / n));
  }
  return s.done();
}
