// Checks of peira/missing.
//
// - Rubin's rules in closed form: estimates 10, 12, 11, 13, 9 with variances 4 each give Q̄ = 11, Ū = 4, B = 2.5,
//   T = 4 + 1.2 × 2.5 = 7, γ = 3/7, Rubin's degrees of freedom 4 × 49/9 = 21.78; Barnard and Rubin's tend to Rubin's as
//   the complete data's grow, and are smaller for small complete data. With identical imputations, B = 0 and T = Ū.
// - Multiple imputation is valid under MAR: 1,500 simulated samples of 400 units, Y depending on a class X that is always
//   seen, and Y missing more often in one class (MAR given X). The complete-case mean is biased; twenty imputations by the
//   approximate Bayesian bootstrap within classes, pooled by Rubin's rules, cover the mean near 95% of the time.

#include <peira/missing/missing.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace peira::missing;

int main() {
  calx::verification::Suite s;

  {
    const std::vector<double> q{10, 12, 11, 13, 9}, u(5, 4.0);
    const Pooled p = pool(q, u);
    s.check("rubin: Q̄ = 11", p.estimate - 11, 1e-14);
    s.check("rubin: B = 2.5", p.between - 2.5, 1e-14);
    s.check("rubin: T = Ū + (1 + 1/m)B = 7", p.total - 7, 1e-14);
    s.check("rubin: γ = 3/7", p.fraction_missing - 3.0 / 7, 1e-15);
    s.check("rubin: Rubin's degrees of freedom (m − 1)/γ² = 21.78", p.df - 4 * 49.0 / 9, 1e-12);
    const Pooled big = pool(q, u, 1e12), small = pool(q, u, 20);
    s.check("barnard-rubin: tends to Rubin's as the complete data's degrees of freedom grow (10¹²)", big.df - p.df, 1e-8);
    s.require("barnard-rubin: smaller with 20 complete-data degrees of freedom", small.df < p.df && small.df < 20, std::to_string(small.df));
    const std::vector<double> same{5, 5, 5};
    const Pooled z = pool(same, std::vector<double>{2, 2, 2});
    s.require("rubin: identical imputations give B = 0 and T = Ū", z.between == 0 && z.total == 2);
  }
  {
    calx::random::Stream st(53, 9, 0, 0);
    const int reps = 1500, n = 400, m = 20;
    int covered = 0;
    double cc_bias = 0;
    // Classes X ∈ {0, 1} with probability 1/2; Y | X ~ N(10 + 6X, 2²); the population mean 13. Missing with probability
    // 0.1 when X = 0 and 0.6 when X = 1 (MAR given X).
    for (int r = 0; r < reps; r++) {
      std::vector<double> y(n);
      std::vector<int> x(n), seen(n);
      double cc = 0, ncc = 0;
      for (int i = 0; i < n; i++) {
        x[i] = static_cast<int>(st.below(2));
        y[i] = st.normal(10 + 6 * x[i], 2);
        seen[i] = st.uniform() >= (x[i] ? 0.6 : 0.1);
        if (seen[i]) { cc += y[i]; ncc++; }
      }
      cc_bias += (cc / ncc - 13) / reps;
      std::vector<double> est(m), var(m);
      for (int l = 0; l < m; l++) {
        const std::vector<double> yi = approximate_bayesian_bootstrap(y, seen, x, st);
        double mean = 0;
        for (double v : yi) mean += v / n;
        double ss = 0;
        for (double v : yi) ss += (v - mean) * (v - mean);
        est[l] = mean, var[l] = ss / (n - 1) / n;
      }
      const Pooled p = pool(est, var, n - 1);
      covered += p.ci_lower <= 13 && 13 <= p.ci_upper;
    }
    const double share = static_cast<double>(covered) / reps;
    s.require("complete cases: biased under MAR (the mean low by about 1)", cc_bias < -0.8, std::to_string(cc_bias));
    s.require("multiple imputation: covers the mean within four binomial standard errors of 95%", std::fabs(share - 0.95) <= 4 * std::sqrt(0.0475 / reps), std::to_string(share));
  }
  return s.done();
}
