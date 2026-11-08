// Checks of peira/design: assignment mechanisms and power.
//
// - Complete randomisation treats exactly n₁ units, and every one of the C(5, 2) = 10 subsets is equally likely (a χ² of
//   20,000 draws on 9 degrees of freedom, below its 0.999 quantile 27.88); each_complete visits the C(10, 4) = 210
//   assignments, each once.
// - Blocked designs treat exactly the number asked in each block; clustered designs keep a cluster in one arm.
// - Bernoulli enumeration: the probabilities of the 2ⁿ assignments sum to one, and E[Z_i] = e_i.
// - Power: δ = σ/2 at α = 0.05 and power 0.8 needs 2 (1.959964 + 0.841621)² / 0.25 = 62.79, so 63 per arm, and 63
//   reaches 0.8 where 62 does not; the power of 30 per arm at δ = 0.6σ agrees with 40,000 simulated z-tests within four
//   standard errors. Proportions 0.6 against 0.4 need 96.92, so 97 per arm (Fleiss et al., eq. 4.14).

#include <peira/design/design.hpp>
#include <peira/design/power.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace peira::design;

int main() {
  calx::verification::Suite s;

  {
    Stream st(7, 1, 0, 0);
    std::map<int, int> counts;
    bool sizes = true;
    const int draws = 20000;
    for (int r = 0; r < draws; r++) {
      const Assignment z = complete(5, 2, st);
      int key = 0, treated = 0;
      for (int i = 0; i < 5; i++) { key = key * 2 + z[i]; treated += z[i]; }
      sizes = sizes && treated == 2;
      counts[key]++;
    }
    double chi2 = 0;
    for (const auto& [k, c] : counts) chi2 += (c - draws / 10.0) * (c - draws / 10.0) / (draws / 10.0);
    s.require("complete: exactly n1 of n treated in every draw", sizes);
    s.require("complete: all ten subsets of size 2 of 5 drawn", counts.size() == 10, std::to_string(counts.size()));
    s.check("complete: subsets equally likely (χ² on 9 df, below its 0.999 quantile 27.88)", chi2, 27.88);
  }
  {
    std::set<std::vector<int>> seen;
    std::size_t visits = 0;
    bool sizes = true;
    each_complete(10, 4, [&](const Assignment& z) {
      visits++;
      seen.insert(z);
      int t = 0;
      for (int zi : z) t += zi;
      sizes = sizes && t == 4;
    });
    s.require("each_complete: C(10, 4) = 210 assignments, each once, each with 4 treated", visits == 210 && seen.size() == 210 && sizes, std::to_string(visits));
    s.check("binomial_coefficient: C(30, 15) = 155117520", binomial_coefficient(30, 15) - 155117520.0, 0);
  }
  {
    Stream st(11, 2, 0, 0);
    const std::vector<int> block{0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2};
    bool ok = true;
    for (int r = 0; r < 200; r++) {
      const Assignment z = blocked(block, {{0, 2}, {1, 1}, {2, 3}}, st);
      int t[3] = {0, 0, 0};
      for (std::size_t i = 0; i < z.size(); i++) t[block[i]] += z[i];
      ok = ok && t[0] == 2 && t[1] == 1 && t[2] == 3;
    }
    s.require("blocked: the number asked treated in each block, every time", ok);
    const std::vector<int> cluster{5, 5, 7, 7, 7, 9, 9, 11, 11, 11};
    bool whole = true, count = true;
    for (int r = 0; r < 200; r++) {
      const Assignment z = clustered(cluster, 2, st);
      std::map<int, std::set<int>> arms;
      for (std::size_t i = 0; i < z.size(); i++) arms[cluster[i]].insert(z[i]);
      int treated = 0;
      for (const auto& [c, a] : arms) { whole = whole && a.size() == 1; treated += *a.begin(); }
      count = count && treated == 2;
    }
    s.require("clustered: every cluster in one arm, two clusters treated", whole && count);
  }
  {
    const std::vector<double> e{0.1, 0.5, 0.8, 0.3, 0.65};
    double total = 0;
    std::vector<double> mean(5, 0.0);
    each_bernoulli(e, [&](const Assignment& z, double p) {
      total += p;
      for (int i = 0; i < 5; i++) mean[i] += p * z[i];
    });
    double worst = 0;
    for (int i = 0; i < 5; i++) worst = std::max(worst, std::fabs(mean[i] - e[i]));
    s.check("each_bernoulli: the probabilities of the 32 assignments sum to one", total - 1, 1e-14);
    s.check("each_bernoulli: E[Z_i] = e_i", worst, 1e-14);
    Stream st(13, 3, 0, 0);
    double treated = 0;
    for (int r = 0; r < 100000; r++) treated += bernoulli(1, 0.3, st)[0];
    s.check("bernoulli: the share treated is p (0.3), within four standard errors", std::fabs(treated / 1e5 - 0.3) / std::sqrt(0.21 / 1e5), 4);
  }
  {
    const std::size_t n = n_per_arm_means(0.5, 1.0, 0.8, 0.05);
    s.require("power: δ = σ/2 at α = 0.05 and power 0.8 needs 63 per arm (62.79 by the normal approximation)", n == 63, std::to_string(n));
    s.require("power: 63 per arm reaches 0.8, 62 does not", power_means(0.5, 1, 63, 63) >= 0.8 && power_means(0.5, 1, 62, 62) < 0.8);
    // The z-test of 30 against 30 at δ = 0.6σ, simulated.
    Stream st(17, 4, 0, 0);
    const int reps = 40000;
    int rejected = 0;
    for (int r = 0; r < reps; r++) {
      double a = 0, b = 0;
      for (int i = 0; i < 30; i++) { a += st.normal(0.6, 1); b += st.normal(0, 1); }
      const double z = (a / 30 - b / 30) / std::sqrt(2.0 / 30);
      rejected += std::fabs(z) > 1.959963984540054;
    }
    const double p = power_means(0.6, 1, 30, 30), sim = static_cast<double>(rejected) / reps;
    s.check("power: the formula against 40,000 simulated z-tests (standard errors)", std::fabs(sim - p) / std::sqrt(p * (1 - p) / reps), 4);
    const std::size_t np = n_per_arm_proportions(0.6, 0.4, 0.8, 0.05);
    s.require("power: proportions 0.6 against 0.4 need 97 per arm (96.92)", np == 97, std::to_string(np));
    s.require("power: 97 per arm reaches 0.8 for proportions, 96 does not", power_proportions(0.6, 0.4, 97) >= 0.8 && power_proportions(0.6, 0.4, 96) < 0.8);
    s.check("design effect: clusters of 20 at ρ = 0.05 multiply the variance by 1.95", design_effect(20, 0.05) - 1.95, 1e-15);
  }
  return s.done();
}
