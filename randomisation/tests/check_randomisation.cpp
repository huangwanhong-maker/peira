// Checks of peira/randomisation, against truths computed by enumerating the design itself.
//
// - A population of eight units with known potential outcomes and a heterogeneous effect, under complete randomisation of
//   four: over all 70 assignments, E[τ̂] = τ exactly; Var(τ̂) = S₁²/4 + S₀²/4 − S_τ²/8 exactly (Neyman 1923); and
//   E[V̂] = S₁²/4 + S₀²/4, strictly above Var(τ̂). With a constant effect the two agree exactly.
// - Blocks: the stratified estimator is unbiased over the product of the blocks' designs (6 × 10 = 60 assignments).
// - Fisher's lady tasting tea (Fisher 1935, §2): eight cups, four with the milk first, all named rightly: the one-sided
//   p-value is 1/C(8, 4) = 1/70.
// - The exact test is of level α: over every assignment of a population where the sharp null holds, the share with
//   p ≤ 0.05 is at most 0.05. Drawn p-values agree with enumerated ones within their Monte Carlo error.
// - Neyman's intervals cover the average effect near 95% of the time over 4,000 simulated experiments of 200 units.

#include <peira/randomisation/randomisation.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace peira::randomisation;
using peira::design::Assignment;

int main() {
  calx::verification::Suite s;

  const auto population_check = [&](const std::vector<double>& y1, const std::vector<double>& y0, const std::string& name, bool constant) {
    const std::size_t n = y1.size();
    double tau = 0, m1 = 0, m0 = 0;
    for (std::size_t i = 0; i < n; i++) { tau += (y1[i] - y0[i]) / n; m1 += y1[i] / n; m0 += y0[i] / n; }
    double S1 = 0, S0 = 0, St = 0;
    for (std::size_t i = 0; i < n; i++) {
      S1 += (y1[i] - m1) * (y1[i] - m1) / (n - 1);
      S0 += (y0[i] - m0) * (y0[i] - m0) / (n - 1);
      St += (y1[i] - y0[i] - tau) * (y1[i] - y0[i] - tau) / (n - 1);
    }
    double mean = 0, sq = 0, ev = 0, count = 0;
    peira::design::each_complete(n, n / 2, [&](const Assignment& z) {
      std::vector<double> y(n);
      for (std::size_t i = 0; i < n; i++) y[i] = z[i] ? y1[i] : y0[i];
      const Neyman r = difference_in_means(y, z);
      mean += r.estimate, sq += r.estimate * r.estimate, ev += r.variance, count++;
    });
    mean /= count, sq /= count, ev /= count;
    const double var = sq - mean * mean, n1 = n / 2.0;
    s.check(name + ": E[τ̂] = τ over every assignment", mean - tau, 1e-12);
    s.check(name + ": Var(τ̂) = S₁²/n₁ + S₀²/n₀ − S_τ²/n exactly", var - (S1 / n1 + S0 / n1 - St / n), 1e-12);
    s.check(name + ": E[V̂] = S₁²/n₁ + S₀²/n₀", ev - (S1 / n1 + S0 / n1), 1e-12);
    if (constant) s.check(name + ": with a constant effect, E[V̂] = Var(τ̂)", ev - var, 1e-12);
    else s.require(name + ": with a heterogeneous effect, E[V̂] > Var(τ̂) (conservative)", ev > var + 1e-6, std::to_string(ev) + " > " + std::to_string(var));
  };
  population_check({12, 9, 15, 7, 11, 14, 6, 10}, {8, 9, 10, 3, 11, 8, 5, 4}, "Neyman, heterogeneous effects", false);
  population_check({11, 10, 13, 6, 14, 10, 8, 7}, {8, 7, 10, 3, 11, 7, 5, 4}, "Neyman, a constant effect", true);

  {
    // Two blocks of four and five units, two treated in each.
    const std::vector<double> y1{5, 7, 6, 9, 20, 24, 21, 26, 22}, y0{3, 6, 4, 5, 15, 18, 19, 20, 17};
    const std::vector<int> block{0, 0, 0, 0, 1, 1, 1, 1, 1};
    double tau = 0;
    for (int i = 0; i < 9; i++) tau += (y1[i] - y0[i]) / 9;
    double mean = 0, count = 0;
    peira::design::each_complete(4, 2, [&](const Assignment& za) {
      peira::design::each_complete(5, 2, [&](const Assignment& zb) {
        std::vector<int> z(za.begin(), za.end());
        z.insert(z.end(), zb.begin(), zb.end());
        std::vector<double> y(9);
        for (int i = 0; i < 9; i++) y[i] = z[i] ? y1[i] : y0[i];
        mean += stratified(y, z, block).estimate, count++;
      });
    });
    s.check("stratified: unbiased over the 60 assignments of two blocks", mean / count - tau, 1e-12);
  }

  {
    // The lady tasting tea: z = the cups with milk first, y = the cups she names so.
    const std::vector<int> z{1, 1, 1, 1, 0, 0, 0, 0};
    const std::vector<double> y{1, 1, 1, 1, 0, 0, 0, 0};
    const FisherTest t = fisher(y, z, 0, Alternative::greater);
    s.check("fisher: the lady tasting tea, all eight right, p = 1/70", t.p_value - 1.0 / 70, 1e-15);
    s.require("fisher: enumerated, over the 70 assignments", t.exact && t.assignments == 70);
  }
  {
    // Level: a population where Y(1) = Y(0) + 2 for every unit; the test of τ₀ = 2 over every assignment.
    const std::vector<double> y0{3.1, 4.7, 2.2, 5.9, 4.4, 3.8, 6.1, 2.9, 5.0, 4.1};
    std::size_t rejected = 0, total = 0;
    peira::design::each_complete(10, 5, [&](const Assignment& z) {
      std::vector<double> y(10);
      for (int i = 0; i < 10; i++) y[i] = y0[i] + 2 * z[i];
      rejected += fisher(y, z, 2).p_value <= 0.05;
      total++;
    });
    s.require("fisher: of level 0.05 when its sharp null holds, over all 252 assignments", static_cast<double>(rejected) / total <= 0.05,
              std::to_string(rejected) + " of " + std::to_string(total));
    // Drawn against enumerated.
    const std::vector<int> z{1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1};
    const std::vector<double> y{6.2, 4.1, 5.9, 7.3, 4.4, 5.0, 6.8, 3.9, 5.5, 4.8, 6.1, 5.2, 4.0, 6.6};
    const FisherTest exact = fisher(y, z, 0), drawn = fisher(y, z, 0, Alternative::two_sided, 0, 40000, 99);
    s.require("fisher: C(14, 7) = 3432 enumerated; drawn when asked", exact.exact && exact.assignments == 3432 && !drawn.exact);
    s.check("fisher: drawn p against the exact (Monte Carlo standard errors)", std::fabs(drawn.p_value - exact.p_value) / std::sqrt(exact.p_value * (1 - exact.p_value) / 40000 + 1e-12), 4);
  }
  {
    // Coverage of Neyman's interval over simulated experiments.
    calx::random::Stream st(23, 5, 0, 0);
    const int reps = 4000, n = 200;
    int covered = 0;
    for (int r = 0; r < reps; r++) {
      std::vector<double> y1(n), y0(n);
      double tau = 0;
      for (int i = 0; i < n; i++) { y0[i] = st.normal(10, 3); y1[i] = y0[i] + st.normal(1.5, 1.0); tau += (y1[i] - y0[i]) / n; }
      const Assignment z = peira::design::complete(n, n / 2, st);
      std::vector<double> y(n);
      for (int i = 0; i < n; i++) y[i] = z[i] ? y1[i] : y0[i];
      const Neyman e = difference_in_means(y, z);
      covered += e.ci_lower <= tau && tau <= e.ci_upper;
    }
    const double share = static_cast<double>(covered) / reps;
    s.require("neyman: the 95% interval covers τ at least 95% of the time less four standard errors (conservative)", share >= 0.95 - 4 * std::sqrt(0.0475 / reps), std::to_string(share));
  }
  return s.done();
}
