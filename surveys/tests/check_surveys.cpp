// Checks of peira/surveys, by enumerating every sample a design can draw.
//
// - Stratified simple random sampling without replacement, two strata of 4 and 5 units, 2 drawn from each (6 × 10 = 60
//   samples, equally likely): the total is unbiased; its variance is Σ N_h²(1 − f_h)S_h²/n_h exactly; and the
//   ultimate-cluster estimator's expectation is that variance, exactly.
// - PSUs drawn with replacement with unequal probabilities (four PSUs, three draws, all 4³ ordered samples): the
//   Hansen–Hurwitz estimator is unbiased, and the ultimate-cluster estimator with f = 0 is unbiased for its variance.
// - One stratum of n units without a finite population correction: the mean's standard error is s/√n.
// - Raking on two margins meets both totals; raking on one margin is post-stratification, in one step.

#include <peira/surveys/surveys.hpp>

#include <peira/design/design.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace peira::surveys;

int main() {
  calx::verification::Suite s;

  {
    const std::vector<double> A{3, 8, 5, 12}, B{20, 14, 25, 17, 22};
    double population = 0;
    for (double v : A) population += v;
    for (double v : B) population += v;
    const auto S2 = [](const std::vector<double>& v) {
      double m = 0;
      for (double x : v) m += x / v.size();
      double ss = 0;
      for (double x : v) ss += (x - m) * (x - m);
      return ss / (v.size() - 1);
    };
    const double var_formula = 16 * (1 - 0.5) * S2(A) / 2 + 25 * (1 - 0.4) * S2(B) / 2;
    double mean = 0, sq = 0, ev = 0, count = 0;
    peira::design::each_complete(4, 2, [&](const peira::design::Assignment& za) {
      peira::design::each_complete(5, 2, [&](const peira::design::Assignment& zb) {
        std::vector<double> y;
        std::vector<int> h;
        for (int i = 0; i < 4; i++) if (za[i]) { y.push_back(A[i]); h.push_back(0); }
        for (int i = 0; i < 5; i++) if (zb[i]) { y.push_back(B[i]); h.push_back(1); }
        const Design d = stratified_srswor(h, {{0, 4}, {1, 5}});
        const Estimate e = total(y, d);
        mean += e.value, sq += e.value * e.value, ev += e.se * e.se, count++;
      });
    });
    mean /= count, sq /= count, ev /= count;
    s.check("stratified SRSWOR: the total is unbiased over all 60 samples", mean - population, 1e-10);
    s.check("stratified SRSWOR: Var = Σ N_h²(1 − f_h)S_h²/n_h exactly", (sq - mean * mean) - var_formula, 1e-9);
    s.check("stratified SRSWOR: the ultimate-cluster estimator is unbiased for it", ev - var_formula, 1e-9);
  }

  {
    // PSUs with totals T and draw probabilities p; three draws with replacement; the units of a PSU drawn carry weight
    // 1/(m p) (here a PSU's total is carried by one unit per draw, the draw its PSU label).
    const std::vector<double> T{40, 25, 70, 15}, p{0.3, 0.2, 0.4, 0.1};
    const double population = 150;
    double mean = 0, sq = 0, ev = 0;
    for (int a = 0; a < 4; a++)
      for (int b = 0; b < 4; b++)
        for (int c = 0; c < 4; c++) {
          const int draw[3] = {a, b, c};
          const double prob = p[a] * p[b] * p[c];
          Design d;
          std::vector<double> y;
          for (int j = 0; j < 3; j++) {
            y.push_back(T[draw[j]]);
            d.weight.push_back(1 / (3 * p[draw[j]]));
            d.psu.push_back(j);
          }
          const Estimate e = peira::surveys::total(y, d);
          mean += prob * e.value, sq += prob * e.value * e.value, ev += prob * e.se * e.se;
        }
    s.check("hansen-hurwitz: unbiased over all 64 ordered samples", mean - population, 1e-10);
    s.check("hansen-hurwitz: the ultimate-cluster estimator (f = 0) is unbiased for its variance", ev - (sq - mean * mean), 1e-9);
  }

  {
    const std::vector<double> y{4.2, 5.1, 3.3, 6.0, 5.5, 4.8, 3.9};
    Design d;
    d.weight.assign(7, 10.0);
    const Estimate e = mean(y, d);
    double m = 0;
    for (double v : y) m += v / 7;
    double ss = 0;
    for (double v : y) ss += (v - m) * (v - m);
    s.check("mean: one stratum without fpc, equal weights: the mean", e.value - m, 1e-14);
    s.check("mean: its standard error is s/√n", e.se - std::sqrt(ss / 6 / 7), 1e-14);
    s.require("mean: on n − 1 = 6 degrees of freedom", e.df == 6);
  }

  {
    const std::vector<double> w(8, 1.0);
    const std::vector<int> sex{0, 0, 0, 1, 1, 1, 1, 0}, age{0, 1, 2, 0, 1, 2, 2, 1};
    const Raking r = rake(w, {sex, age}, {{{0, 48}, {1, 52}}, {{0, 30}, {1, 40}, {2, 30}}});
    s.require("rake: converges on two margins", r.converged, std::to_string(r.iterations) + " iterations, gap " + std::to_string(r.gap));
    double female = 0, young = 0;
    for (int i = 0; i < 8; i++) { female += sex[i] ? r.weight[i] : 0; young += age[i] == 0 ? r.weight[i] : 0; }
    s.check("rake: the first margin met", female - 52, 1e-8);
    s.check("rake: the second margin met", young - 30, 1e-8);
    const Raking post = rake(w, {sex}, {{{0, 48}, {1, 52}}});
    s.require("rake: one margin is post-stratification, reached in one step", post.converged && post.iterations == 1);
    s.check("rake: post-stratified weights are N_g/n_g", post.weight[0] - 12.0, 1e-12);
  }
  return s.done();
}
