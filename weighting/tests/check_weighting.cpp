// Checks of peira/weighting.
//
// - A design of fixed size 2 on five units with unequal probabilities, given as the probability of each of its ten
//   samples: over every sample, E[Σ y/π] is the total exactly, and both variance estimators, Horvitz–Thompson's and
//   Sen–Yates–Grundy's, have the variance of the total as their expectation, exactly.
// - Hájek's mean of a constant is that constant in every sample; Horvitz–Thompson's is not.
// - The inverse-probability-weighted effect (unnormalised) is unbiased for the average effect over all 2⁶ assignments of
//   a Bernoulli design with unit-specific probabilities, exactly.
// - Double robustness, by 300 simulated studies of 2,000 units (x ~ N(0, 1), e(x) = logistic(−0.5 + 1.2x), Y(0) =
//   1 + x + 0.8x² + ε, Y(1) = Y(0) + 2 + 0.5x, so τ = 2): the augmented estimator is unbiased (within four Monte Carlo
//   standard errors) when the propensity is right and the outcome model wrong, and when the outcome model is right and the
//   propensity wrong; when both are wrong it is biased, by more than four. Its intervals cover τ near 95% of the time
//   when both are right.

#include <peira/design/design.hpp>
#include <peira/weighting/weighting.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace peira::weighting;

int main() {
  calx::verification::Suite s;

  {
    const std::vector<double> Y{3, 7, 1, 9, 5};
    double total = 0;
    for (double v : Y) total += v;
    // The design: P({i, j}) ∝ (i + 1) + (j + 1) + 1, over the ten pairs.
    std::vector<std::pair<int, int>> pairs;
    std::vector<double> p;
    double z = 0;
    for (int i = 0; i < 5; i++)
      for (int j = i + 1; j < 5; j++) { pairs.push_back({i, j}); p.push_back(i + j + 3.0); z += i + j + 3.0; }
    for (double& v : p) v /= z;
    std::vector<double> pi(5, 0.0);
    Matrix pij(5);
    for (std::size_t k = 0; k < pairs.size(); k++) {
      const auto [i, j] = pairs[k];
      pi[i] += p[k], pi[j] += p[k];
      pij(i, j) += p[k], pij(j, i) += p[k];
    }
    double mean = 0, sq = 0, eht = 0, esyg = 0, hajek_const = 0, ht_const_spread = 0;
    for (std::size_t k = 0; k < pairs.size(); k++) {
      const auto [i, j] = pairs[k];
      const std::vector<double> y{Y[i], Y[j]}, pis{pi[i], pi[j]};
      Matrix ps(2);
      ps(0, 1) = ps(1, 0) = pij(i, j);
      const double t = horvitz_thompson_total(y, pis);
      mean += p[k] * t, sq += p[k] * t * t;
      eht += p[k] * ht_variance(y, pis, ps);
      esyg += p[k] * syg_variance(y, pis, ps);
      const std::vector<double> four{4, 4};
      hajek_const = std::max(hajek_const, std::fabs(hajek_mean(four, pis) - 4));
      ht_const_spread = std::max(ht_const_spread, std::fabs(horvitz_thompson_total(four, pis) / 5 - 4));
    }
    const double var = sq - mean * mean;
    s.check("horvitz-thompson: E[Σ y/π] is the total over the design's ten samples", mean - total, 1e-12);
    s.check("horvitz-thompson: E[V̂_HT] = Var exactly", eht - var, 1e-10);
    s.check("sen-yates-grundy: E[V̂_SYG] = Var exactly (a design of fixed size)", esyg - var, 1e-10);
    s.check("hajek: the mean of a constant is that constant in every sample", hajek_const, 1e-14);
    s.require("horvitz-thompson: the mean of a constant is not (the weights' sum varies)", ht_const_spread > 0.01, std::to_string(ht_const_spread));
  }

  {
    const std::vector<double> e{0.2, 0.5, 0.7, 0.35, 0.6, 0.45}, y1{5, 8, 6, 9, 4, 7}, y0{3, 6, 6, 4, 2, 5};
    double tau = 0;
    for (int i = 0; i < 6; i++) tau += (y1[i] - y0[i]) / 6;
    double mean = 0;
    peira::design::each_bernoulli(e, [&](const peira::design::Assignment& zz, double p) {
      std::vector<double> y(6);
      for (int i = 0; i < 6; i++) y[i] = zz[i] ? y1[i] : y0[i];
      mean += p * ipw(y, zz, e, false).estimate;
    });
    s.check("ipw: unbiased for the average effect over all 64 assignments of a Bernoulli design", mean - tau, 1e-12);
  }

  {
    calx::random::Stream st(31, 6, 0, 0);
    const int reps = 300, n = 2000;
    double sums[3] = {0, 0, 0}, squares[3] = {0, 0, 0};
    int covered = 0;
    for (int r = 0; r < reps; r++) {
      std::vector<double> x(n), y(n);
      std::vector<int> z(n);
      for (int i = 0; i < n; i++) {
        x[i] = st.normal();
        const double e = 1 / (1 + std::exp(-(-0.5 + 1.2 * x[i])));
        z[i] = st.uniform() < e;
        const double y0 = 1 + x[i] + 0.8 * x[i] * x[i] + st.normal();
        y[i] = z[i] ? y0 + 2 + 0.5 * x[i] : y0;
      }
      Rectangular right_p(n, 2), wrong_p(n, 1), right_q(n, 3), wrong_q(n, 2);
      for (int i = 0; i < n; i++) {
        right_p(i, 0) = 1, right_p(i, 1) = x[i];
        wrong_p(i, 0) = 1;
        right_q(i, 0) = 1, right_q(i, 1) = x[i], right_q(i, 2) = x[i] * x[i];
        wrong_q(i, 0) = 1, wrong_q(i, 1) = x[i];
      }
      const Rectangular* P[3] = {&right_p, &wrong_p, &wrong_p};
      const Rectangular* Q[3] = {&wrong_q, &right_q, &wrong_q};
      for (int c = 0; c < 3; c++) {
        const Nuisances nu = fit_nuisances(*P[c], *Q[c], y, z);
        const Effect a = aipw(y, z, nu.e, nu.m1, nu.m0);
        sums[c] += a.estimate, squares[c] += a.estimate * a.estimate;
      }
      const Nuisances both = fit_nuisances(right_p, right_q, y, z);
      const Effect a = aipw(y, z, both.e, both.m1, both.m0);
      covered += a.ci_lower <= 2 && 2 <= a.ci_upper;
    }
    const char* names[3] = {"propensity right, outcome wrong", "outcome right, propensity wrong", "both wrong"};
    for (int c = 0; c < 3; c++) {
      const double m = sums[c] / reps, sd = std::sqrt(squares[c] / reps - m * m), z = std::fabs(m - 2) / (sd / std::sqrt(reps));
      if (c < 2) s.check(std::string("aipw: unbiased when ") + names[c] + " (Monte Carlo standard errors from τ = 2)", z, 4);
      else s.require(std::string("aipw: biased when ") + names[c] + ", by more than four Monte Carlo standard errors", z > 4, "mean " + std::to_string(m) + ", z " + std::to_string(z));
    }
    const double share = static_cast<double>(covered) / reps;
    s.require("aipw: 95% intervals cover τ when both are right, within four binomial standard errors", std::fabs(share - 0.95) <= 4 * std::sqrt(0.0475 / reps), std::to_string(share));
  }
  return s.done();
}
