// Checks of peira/detection: site occupancy with imperfect detection (occupancy.hpp, on latent.hpp).
//
// - The likelihood equals the brute-force sum over every configuration of the latent occupancy (2⁵ of five sites), with
//   covariates and a visit not made.
// - The analytic gradient and the Hessian of Louis's identity agree with central differences.
// - Identification: with one visit, ℓ depends on (ψ, p) only through ψp; a visit not made by any site is the same as a visit
//   that does not exist.
// - The constant model's MLE is the known characterisation: p̂ solves d/(s_D T) = p/(1 − (1 − p)^T) (a binomial truncated at
//   zero, solved here by calx's Brent), and ψ̂ = s_D/(s (1 − (1 − p̂)^T)).
// - The posterior of a site with no detection is ψ(1 − p)^T/(ψ(1 − p)^T + 1 − ψ); the finite-sample occupancy is their mean
//   with the detected sites counted whole.
// - Recovery from 500 simulated surveys (300 sites, 4 visits, a covariate in ψ and one in p): every coefficient's bias below a
//   quarter of the estimates' spread (the maximum likelihood estimate is consistent, with a bias of order 1/s; the Monte Carlo
//   error of this ratio is 0.045), the Wald intervals' coverage within 3.5 Monte Carlo standard errors of 95%, the standard
//   errors' mean within 15% of the estimates' spread; the finite-sample occupancy unbiased for the realised occupancy (within four
//   Monte Carlo standard errors) while the naive one falls short; AIC prefers the covariate that is there.

#include <peira/detection/detection.hpp>

#include <calx/random/stream.hpp>
#include <calx/solvers/roots.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace pd = peira::detection;
using calx::linalg::Rectangular;

namespace {

double logistic(double x) { return 1 / (1 + std::exp(-x)); }

// The likelihood by enumerating z ∈ {0, 1}^S, with plain probabilities.
double brute_force(const pd::SiteVisits& y, const Rectangular& X, const Rectangular& W, const std::vector<double>& theta) {
  const std::size_t S = y.sites, T = y.visits;
  double L = 0;
  for (unsigned z = 0; z < (1u << S); z++) {
    double P = 1;
    for (std::size_t i = 0; i < S; i++) {
      double ex = 0;
      for (std::size_t a = 0; a < X.cols; a++) ex += X(i, a) * theta[a];
      const double psi = logistic(ex);
      const bool zi = (z >> i) & 1u;
      P *= zi ? psi : 1 - psi;
      for (std::size_t j = 0; j < T; j++) {
        const int v = y(i, j);
        if (v < 0) continue;
        double ew = 0;
        for (std::size_t b = 0; b < W.cols; b++) ew += W(i * T + j, b) * theta[X.cols + b];
        const double p = logistic(ew);
        P *= zi ? (v ? p : 1 - p) : (v ? 0.0 : 1.0);
      }
    }
    L += P;
  }
  return std::log(L);
}

double worst_relative(const std::vector<double>& a, const std::vector<double>& b) {
  double w = 0;
  for (std::size_t i = 0; i < a.size(); i++) w = std::max(w, std::fabs(a[i] - b[i]) / std::max(1.0, std::fabs(b[i])));
  return w;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- A tiny survey: five sites, three visits, one not made; ψ on a site covariate, p on a visit covariate.
  pd::SiteVisits tiny(5, 3);
  const int hist[5][3] = {{1, 0, 1}, {0, 0, 0}, {0, 1, -1}, {0, 0, 0}, {1, 1, 1}};
  for (std::size_t i = 0; i < 5; i++)
    for (std::size_t j = 0; j < 3; j++) tiny(i, j) = hist[i][j];
  Rectangular X(5, 2), W(15, 2);
  const double xs[5] = {-1.2, 0.3, 0.8, -0.4, 1.5};
  for (std::size_t i = 0; i < 5; i++) {
    X(i, 0) = 1, X(i, 1) = xs[i];
    for (std::size_t j = 0; j < 3; j++) W(i * 3 + j, 0) = 1, W(i * 3 + j, 1) = 0.5 * static_cast<double>(j) - 0.3 * xs[i];
  }
  const std::vector<double> theta{0.4, -0.7, -0.2, 0.9};
  {
    const double l = pd::occupancy_log_likelihood(tiny, X, W, theta);
    s.check("occupancy: ℓ equals the brute-force sum over the 2⁵ latent occupancies", l - brute_force(tiny, X, W, theta), 1e-13);
    const auto e = pd::detail::latent_evaluate(tiny, X, W, pd::Mixing::bernoulli, 1, theta, 2);
    std::vector<double> num(4);
    for (std::size_t a = 0; a < 4; a++) {
      auto tp = theta, tm = theta;
      const double h = 1e-5;
      tp[a] += h, tm[a] -= h;
      num[a] = (pd::occupancy_log_likelihood(tiny, X, W, tp) - pd::occupancy_log_likelihood(tiny, X, W, tm)) / (2 * h);
    }
    s.check("occupancy: the analytic gradient against central differences", worst_relative(e.gradient, num), 1e-8);
    double worst = 0;
    for (std::size_t a = 0; a < 4; a++) {
      auto tp = theta, tm = theta;
      const double h = 1e-5;
      tp[a] += h, tm[a] -= h;
      const auto gp = pd::detail::latent_evaluate(tiny, X, W, pd::Mixing::bernoulli, 1, tp, 1).gradient;
      const auto gm = pd::detail::latent_evaluate(tiny, X, W, pd::Mixing::bernoulli, 1, tm, 1).gradient;
      for (std::size_t b = 0; b < 4; b++) worst = std::max(worst, std::fabs((gp[b] - gm[b]) / (2 * h) - e.hessian(a, b)));
    }
    s.check("occupancy: the Hessian of Louis's identity against differences of the gradient", worst, 1e-8);
  }
  {
    // One visit: ℓ(ψ, p) = ℓ(ψ', p') whenever ψp = ψ'p'.
    pd::SiteVisits one(6, 1);
    const int v[6] = {1, 0, 0, 1, 0, 1};
    for (std::size_t i = 0; i < 6; i++) one(i, 0) = v[i];
    const auto I = pd::intercept(6);
    auto logit = [](double p) { return std::log(p / (1 - p)); };
    const double a = pd::occupancy_log_likelihood(one, I, I, std::vector<double>{logit(0.8), logit(0.5)});
    const double b = pd::occupancy_log_likelihood(one, I, I, std::vector<double>{logit(0.5), logit(0.8)});
    s.check("identification: with one visit ℓ depends on ψ and p only through ψp", a - b, 1e-13);
    // A visit nobody made is a visit that does not exist.
    pd::SiteVisits with(5, 4), without(5, 3);
    Rectangular W4(20, 2), W3(15, 2);
    for (std::size_t i = 0; i < 5; i++) {
      for (std::size_t j = 0; j < 3; j++) {
        with(i, j) = without(i, j) = hist[i][j];
        for (std::size_t c = 0; c < 2; c++) W4(i * 4 + j, c) = W3(i * 3 + j, c) = W(i * 3 + j, c);
      }
      with(i, 3) = -1;
      W4(i * 4 + 3, 0) = 1, W4(i * 4 + 3, 1) = 7;
    }
    s.check("identification: a visit made by no site changes nothing",
            pd::occupancy_log_likelihood(with, X, W4, theta) - pd::occupancy_log_likelihood(without, X, W3, theta), 1e-14);
  }

  // ---- The constant model, against its characterisation.
  {
    calx::random::Stream st(2002, 1, 0, 0);
    const std::size_t S = 300, T = 4;
    std::vector<double> psi(S, 0.6), p(S * T, 0.35);
    const auto sim = pd::simulate_visits(st, pd::Mixing::bernoulli, psi, p, T);
    const auto fit = pd::occupancy(sim.y, pd::intercept(S), pd::intercept(S * T));
    double sD = 0, d = 0;
    for (std::size_t i = 0; i < S; i++) {
      double yi = 0;
      for (std::size_t j = 0; j < T; j++) yi += sim.y(i, j);
      d += yi, sD += yi > 0;
    }
    const double target = d / (sD * T);
    const auto root = calx::solvers::brent([&](double q) { return q / (1 - std::pow(1 - q, T)) - target; }, 1e-9, 1 - 1e-12);
    const double phat = root.x, psihat = sD / (S * (1 - std::pow(1 - phat, T)));
    s.require("constant occupancy: converges, with a positive definite information", fit.converged && fit.information_positive,
              std::to_string(fit.iterations) + " iterations");
    s.check("constant occupancy: p̂ solves the truncated binomial's equation d/(s_D T) = p/(1 − (1 − p)^T)", logistic(fit.p_coefficients[0]) - phat, 1e-8);
    s.check("constant occupancy: ψ̂ = s_D / (s (1 − (1 − p̂)^T))", logistic(fit.psi_coefficients[0]) - psihat, 1e-8);
    s.check("constant occupancy: the naive occupancy is s_D/s", fit.naive - sD / S, 0);
    // The posterior of a site never detected, and the finite-sample occupancy.
    const double q0 = psihat * std::pow(1 - phat, T), post0 = q0 / (q0 + 1 - psihat);
    double worst = 0, finite = 0;
    for (std::size_t i = 0; i < S; i++) {
      double yi = 0;
      for (std::size_t j = 0; j < T; j++) yi += sim.y(i, j);
      worst = std::max(worst, std::fabs(fit.occupied[i] - (yi > 0 ? 1.0 : post0)));
      finite += yi > 0 ? 1.0 : post0;
    }
    s.check("constant occupancy: P(z = 1 | no detection) = ψ(1 − p)^T / (ψ(1 − p)^T + 1 − ψ)", worst, 1e-8);
    s.check("constant occupancy: the finite-sample occupancy = (s_D + Σ posteriors of the rest)/s", fit.finite_sample - finite / S, 1e-8);
    s.check("constant occupancy: AIC = −2ℓ + 2k", fit.aic - (-2 * fit.log_likelihood + 4), 1e-12);
  }

  // ---- Recovery from simulated surveys with covariates.
  {
    const std::size_t S = 300, T = 4, R = 500;
    const std::vector<double> truth{0.3, 0.8, -0.4, 0.6};
    std::vector<double> sum(4, 0.0), sum2(4, 0.0), se_sum(4, 0.0), cover(4, 0.0);
    double naive = 0, finite = 0, realised = 0, finite_err2 = 0, short2 = 0, aic_wins = 0;
    int failed = 0;
    for (std::size_t rep = 0; rep < R; rep++) {
      calx::random::Stream st(2018, 2, static_cast<std::uint32_t>(rep), 0);
      Rectangular X(S, 2), W(S * T, 2);
      std::vector<double> psi(S), p(S * T);
      for (std::size_t i = 0; i < S; i++) {
        X(i, 0) = 1, X(i, 1) = st.normal();
        psi[i] = logistic(truth[0] + truth[1] * X(i, 1));
        for (std::size_t j = 0; j < T; j++) {
          W(i * T + j, 0) = 1, W(i * T + j, 1) = st.normal();
          p[i * T + j] = logistic(truth[2] + truth[3] * W(i * T + j, 1));
        }
      }
      const auto sim = pd::simulate_visits(st, pd::Mixing::bernoulli, psi, p, T);
      const auto fit = pd::occupancy(sim.y, X, W);
      if (!fit.converged || !fit.information_positive) { failed++; continue; }
      for (std::size_t a = 0; a < 4; a++) {
        sum[a] += fit.theta[a], sum2[a] += fit.theta[a] * fit.theta[a];
        const double se = std::sqrt(fit.covariance(a, a));
        se_sum[a] += se;
        cover[a] += std::fabs(fit.theta[a] - truth[a]) <= 1.959963984540054 * se;
      }
      double z = 0;
      for (int v : sim.latent) z += v;
      naive += fit.naive, finite += fit.finite_sample, realised += z / S;
      finite_err2 += (fit.finite_sample - z / S) * (fit.finite_sample - z / S);
      short2 += (z / S - fit.naive) * (z / S - fit.naive);
      const auto constant = pd::occupancy(sim.y, pd::intercept(S), W);
      aic_wins += fit.aic < constant.aic;
    }
    const double n = static_cast<double>(R - static_cast<std::size_t>(failed));
    s.require("occupancy recovery: every simulated survey converges with a positive definite information", failed == 0, std::to_string(failed) + " failed");
    const char* names[4] = {"ψ's intercept", "ψ's slope", "p's intercept", "p's slope"};
    for (std::size_t a = 0; a < 4; a++) {
      const double mean = sum[a] / n, sd = std::sqrt(std::max(0.0, sum2[a] / n - mean * mean)), mcse = sd / std::sqrt(n);
      s.check(std::string("occupancy recovery: ") + names[a] + "'s bias / the estimates' spread (Monte Carlo error " + calx::verification::shortest(std::round(mcse / sd * 1000) / 1000) + ")",
              (mean - truth[a]) / sd, 0.25);
      s.check(std::string("occupancy recovery: ") + names[a] + "'s Wald coverage − 0.95 (limit 3.5 Monte Carlo s.e.)", cover[a] / n - 0.95,
              3.5 * std::sqrt(0.95 * 0.05 / n));
      s.check(std::string("occupancy recovery: ") + names[a] + "'s mean standard error / the estimates' spread − 1", se_sum[a] / n / sd - 1, 0.15);
    }
    const double gap_mcse = std::sqrt(finite_err2 / n) / std::sqrt(n);
    s.check("occupancy recovery: the finite-sample occupancy unbiased for the realised (in Monte Carlo s.e.)", (finite - realised) / n / gap_mcse, 4);
    const double shortfall = (realised - naive) / n, short_sd = std::sqrt(std::max(0.0, short2 / n - shortfall * shortfall));
    s.require("occupancy recovery: the naive occupancy falls short of the realised by more than 20 Monte Carlo s.e.",
              shortfall > 20 * short_sd / std::sqrt(n), std::to_string(shortfall) + " (s.e. " + std::to_string(short_sd / std::sqrt(n)) + ")");
    s.require("occupancy recovery: AIC prefers the covariate in ψ in at least 95% of surveys", aic_wins / n >= 0.95, std::to_string(aic_wins / n));
  }
  return s.done();
}
