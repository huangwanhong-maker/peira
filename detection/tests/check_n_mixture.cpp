// Checks of peira/detection: abundance from repeated counts, the N-mixture model (n_mixture.hpp, on latent.hpp).
//
// - The likelihood equals the sum over each site's latent count written out with tyche's Poisson and binomial probabilities
//   (and the negative binomial's from calx's log Γ), on a small survey with covariates and a visit not made.
// - The analytic gradient and Louis's Hessian agree with central differences, for the Poisson and the negative binomial.
// - Thinning: with one visit, ℓ = Σ log Poisson(y_i; λp) exactly (the model is not identified), and the posterior of N_i − y_i
//   is Poisson(λ(1 − p)).
// - The negative binomial tends to the Poisson as r → ∞.
// - The bound K: estimates at K and 2K agree when the reported effect of doubling it is negligible; a bound at the largest count
//   is reported as too small and moves the estimates.
// - Recovery from 200 simulated surveys (150 sites, 4 visits, a covariate in λ and one in p): every coefficient's bias below a
//   quarter of the estimates' spread, the Wald coverage within 3.5 Monte Carlo standard errors of 95%, the standard errors'
//   mean within 15% of the spread; the predicted total Σ E[N_i | y_i] unbiased for the realised total.
// - The negative binomial on 1,000 overdispersed sites: r̂ within three standard errors of the truth, and AIC prefers it to the
//   Poisson.

#include <peira/detection/detection.hpp>

#include <tyche/probability/distributions.hpp>

#include <calx/random/stream.hpp>
#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace pd = peira::detection;
using calx::linalg::Rectangular;

namespace {

double logistic(double x) { return 1 / (1 + std::exp(-x)); }

// The likelihood written out site by site with plain probabilities.
double written_out(const pd::SiteVisits& y, const Rectangular& X, const Rectangular& W, int K, pd::Mixing m, const std::vector<double>& theta) {
  const std::size_t T = y.visits;
  double l = 0;
  for (std::size_t i = 0; i < y.sites; i++) {
    double ex = 0;
    for (std::size_t a = 0; a < X.cols; a++) ex += X(i, a) * theta[a];
    const double lambda = std::exp(ex);
    double Li = 0;
    for (int N = 0; N <= K; N++) {
      double P;
      if (m == pd::Mixing::poisson) {
        P = tyche::probability::Poisson{lambda}.pmf(N);
      } else {
        const double r = std::exp(theta.back());
        P = std::exp(calx::special_functions::log_gamma(N + r) - calx::special_functions::log_gamma(r) - calx::special_functions::log_gamma(N + 1.0) +
                     r * std::log(r / (r + lambda)) + N * std::log(lambda / (r + lambda)));
      }
      for (std::size_t j = 0; j < T; j++) {
        if (y(i, j) < 0) continue;
        double ew = 0;
        for (std::size_t b = 0; b < W.cols; b++) ew += W(i * T + j, b) * theta[X.cols + b];
        P *= tyche::probability::Binomial{N, logistic(ew)}.pmf(y(i, j));
      }
      Li += P;
    }
    l += std::log(Li);
  }
  return l;
}

double derivative_errors(const pd::SiteVisits& y, const Rectangular& X, const Rectangular& W, int K, pd::Mixing m, const std::vector<double>& theta,
                         double& hessian_error) {
  const auto e = pd::detail::latent_evaluate(y, X, W, m, K, theta, 2);
  double gw = 0, hw = 0;
  for (std::size_t a = 0; a < theta.size(); a++) {
    auto tp = theta, tm = theta;
    const double h = 1e-5;
    tp[a] += h, tm[a] -= h;
    const double num = (pd::n_mixture_log_likelihood(y, X, W, K, m, tp) - pd::n_mixture_log_likelihood(y, X, W, K, m, tm)) / (2 * h);
    gw = std::max(gw, std::fabs(num - e.gradient[a]) / std::max(1.0, std::fabs(num)));
    const auto gp = pd::detail::latent_evaluate(y, X, W, m, K, tp, 1).gradient, gm = pd::detail::latent_evaluate(y, X, W, m, K, tm, 1).gradient;
    for (std::size_t b = 0; b < theta.size(); b++) hw = std::max(hw, std::fabs((gp[b] - gm[b]) / (2 * h) - e.hessian(a, b)) / std::max(1.0, std::fabs(e.hessian(a, b))));
  }
  hessian_error = hw;
  return gw;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- A small survey: four sites, three visits (one not made), λ on a site covariate, p on a visit covariate.
  pd::SiteVisits small(4, 3);
  const int counts[4][3] = {{3, 1, 2}, {0, 0, 1}, {5, -1, 4}, {0, 0, 0}};
  for (std::size_t i = 0; i < 4; i++)
    for (std::size_t j = 0; j < 3; j++) small(i, j) = counts[i][j];
  Rectangular X(4, 2), W(12, 2);
  const double xs[4] = {0.5, -1.0, 1.3, -0.2};
  for (std::size_t i = 0; i < 4; i++) {
    X(i, 0) = 1, X(i, 1) = xs[i];
    for (std::size_t j = 0; j < 3; j++) W(i * 3 + j, 0) = 1, W(i * 3 + j, 1) = 0.4 * static_cast<double>(j) - 0.2;
  }
  const std::vector<double> tp{1.1, 0.4, 0.3, -0.6}, tn{1.1, 0.4, 0.3, -0.6, std::log(2.5)};
  const int K = 60;
  s.check("n-mixture (Poisson): ℓ equals the sum over N written out with tyche's probabilities",
          pd::n_mixture_log_likelihood(small, X, W, K, pd::Mixing::poisson, tp) - written_out(small, X, W, K, pd::Mixing::poisson, tp), 1e-12);
  s.check("n-mixture (negative binomial): ℓ equals the sum over N written out",
          pd::n_mixture_log_likelihood(small, X, W, K, pd::Mixing::negative_binomial, tn) - written_out(small, X, W, K, pd::Mixing::negative_binomial, tn), 1e-11);
  {
    double he = 0;
    const double ge = derivative_errors(small, X, W, K, pd::Mixing::poisson, tp, he);
    s.check("n-mixture (Poisson): the analytic gradient against central differences", ge, 1e-7);
    s.check("n-mixture (Poisson): Louis's Hessian against differences of the gradient", he, 1e-7);
    const double gn = derivative_errors(small, X, W, K, pd::Mixing::negative_binomial, tn, he);
    s.check("n-mixture (negative binomial): the analytic gradient against central differences", gn, 1e-7);
    s.check("n-mixture (negative binomial): Louis's Hessian against differences of the gradient", he, 1e-7);
  }
  {
    // r → ∞: the negative binomial's likelihood tends to the Poisson's.
    auto big = tn;
    big.back() = std::log(1e9);
    s.check("n-mixture: the negative binomial with r = 10⁹ gives the Poisson likelihood",
            pd::n_mixture_log_likelihood(small, X, W, K, pd::Mixing::negative_binomial, big) - pd::n_mixture_log_likelihood(small, X, W, K, pd::Mixing::poisson, tp), 1e-6);
  }
  {
    // One visit: thinning.
    calx::random::Stream st(2004, 1, 0, 0);
    const std::size_t S = 40;
    std::vector<double> lam(S, 6.0), p(S, 0.3);
    const auto sim = pd::simulate_visits(st, pd::Mixing::poisson, lam, p, 1);
    const auto I = pd::intercept(S);
    const std::vector<double> th{std::log(6.0), std::log(0.3 / 0.7)};
    double direct = 0;
    for (std::size_t i = 0; i < S; i++) direct += tyche::probability::Poisson{6.0 * 0.3}.log_pmf(sim.y(i, 0));
    s.check("thinning: with one visit, ℓ = Σ log Poisson(y; λp)", pd::n_mixture_log_likelihood(sim.y, I, I, 80, pd::Mixing::poisson, th) - direct, 1e-10);
    std::vector<std::vector<double>> post;
    pd::detail::latent_evaluate(sim.y, I, I, pd::Mixing::poisson, 80, th, 0, &post);
    double worst = 0;
    for (std::size_t i = 0; i < S; i++)
      for (int N = sim.y(i, 0); N <= 80; N++) worst = std::max(worst, std::fabs(post[i][N] - tyche::probability::Poisson{6.0 * 0.7}.pmf(N - sim.y(i, 0))));
    s.check("thinning: the posterior of N − y is Poisson(λ(1 − p))", worst, 1e-12);
  }

  // ---- The bound K.
  {
    calx::random::Stream st(2004, 2, 0, 0);
    const std::size_t S = 120, T = 3;
    std::vector<double> lam(S, 5.0), p(S * T, 0.4);
    const auto sim = pd::simulate_visits(st, pd::Mixing::poisson, lam, p, T);
    int ymax = 0;
    for (int v : sim.y.y) ymax = std::max(ymax, v);
    const auto I = pd::intercept(S), J = pd::intercept(S * T);
    const auto a = pd::n_mixture(sim.y, I, J, 50), b = pd::n_mixture(sim.y, I, J, 100), c = pd::n_mixture(sim.y, I, J, ymax);
    s.require("the bound K: converges at K = 50", a.converged && a.information_positive, std::to_string(a.iterations) + " iterations");
    s.check("the bound K: doubling K = 50 changes ℓ at the estimate by (reported)", a.truncation_effect, 1e-10);
    s.check("the bound K: the prior mass beyond K = 50 (reported)", a.tail_probability, 1e-12);
    s.check("the bound K: λ̂ at K = 50 and K = 100 agree", a.lambda_coefficients[0] - b.lambda_coefficients[0], 1e-8);
    s.check("the bound K: p̂ at K = 50 and K = 100 agree", a.p_coefficients[0] - b.p_coefficients[0], 1e-8);
    s.require("the bound K: at the largest count (" + std::to_string(ymax) + ") it is reported too small and moves λ̂",
              c.truncation_effect > 1e-3 && c.tail_probability > 1e-3 && std::fabs(c.lambda_coefficients[0] - a.lambda_coefficients[0]) > 1e-3,
              "effect " + calx::verification::exponential(c.truncation_effect) + ", tail " + calx::verification::exponential(c.tail_probability));
    const auto pr = pd::predict_lambda(a, std::vector<double>{1.0});
    s.check("predict_lambda: exp(α̂), with the delta method's standard error λ̂ se(α̂)",
            std::fabs(pr.estimate - std::exp(a.lambda_coefficients[0])) + std::fabs(pr.se - pr.estimate * a.lambda_se[0]), 1e-12);
  }

  // ---- Recovery from simulated surveys with covariates (Poisson).
  {
    const std::size_t S = 150, T = 4, R = 200;
    const std::vector<double> truth{std::log(4.0), 0.5, 0.2, -0.5};
    std::vector<double> sum(4, 0.0), sum2(4, 0.0), se_sum(4, 0.0), cover(4, 0.0);
    double gap = 0, gap2 = 0, worst_effect = 0;
    int failed = 0;
    for (std::size_t rep = 0; rep < R; rep++) {
      calx::random::Stream st(2004, 3, static_cast<std::uint32_t>(rep), 0);
      Rectangular Xr(S, 2), Wr(S * T, 2);
      std::vector<double> lam(S), p(S * T);
      for (std::size_t i = 0; i < S; i++) {
        Xr(i, 0) = 1, Xr(i, 1) = st.normal();
        lam[i] = std::exp(truth[0] + truth[1] * Xr(i, 1));
        for (std::size_t j = 0; j < T; j++) {
          Wr(i * T + j, 0) = 1, Wr(i * T + j, 1) = st.normal();
          p[i * T + j] = logistic(truth[2] + truth[3] * Wr(i * T + j, 1));
        }
      }
      const auto sim = pd::simulate_visits(st, pd::Mixing::poisson, lam, p, T);
      const auto fit = pd::n_mixture(sim.y, Xr, Wr, 80);
      if (!fit.converged || !fit.information_positive) { failed++; continue; }
      worst_effect = std::max(worst_effect, fit.truncation_effect);
      for (std::size_t a = 0; a < 4; a++) {
        sum[a] += fit.theta[a], sum2[a] += fit.theta[a] * fit.theta[a];
        const double se = std::sqrt(fit.covariance(a, a));
        se_sum[a] += se;
        cover[a] += std::fabs(fit.theta[a] - truth[a]) <= 1.959963984540054 * se;
      }
      double total = 0;
      for (int v : sim.latent) total += v;
      gap += fit.total_predicted - total, gap2 += (fit.total_predicted - total) * (fit.total_predicted - total);
    }
    const double n = static_cast<double>(R - static_cast<std::size_t>(failed));
    s.require("n-mixture recovery: every simulated survey converges with a positive definite information", failed == 0, std::to_string(failed) + " failed");
    s.check("n-mixture recovery: the largest effect of doubling K = 80 (reported)", worst_effect, 1e-8);
    const char* names[4] = {"λ's intercept", "λ's slope", "p's intercept", "p's slope"};
    for (std::size_t a = 0; a < 4; a++) {
      const double mean = sum[a] / n, sd = std::sqrt(std::max(0.0, sum2[a] / n - mean * mean));
      s.check(std::string("n-mixture recovery: ") + names[a] + "'s bias / the estimates' spread (Monte Carlo error 0.071)", (mean - truth[a]) / sd, 0.3);
      s.check(std::string("n-mixture recovery: ") + names[a] + "'s Wald coverage − 0.95 (limit 3.5 Monte Carlo s.e.)", cover[a] / n - 0.95, 3.5 * std::sqrt(0.95 * 0.05 / n));
      s.check(std::string("n-mixture recovery: ") + names[a] + "'s mean standard error / the estimates' spread − 1", se_sum[a] / n / sd - 1, 0.15);
    }
    const double m = gap / n, sd = std::sqrt(std::max(0.0, gap2 / n - m * m));
    s.check("n-mixture recovery: Σ E[N_i | y_i] unbiased for the realised total (in Monte Carlo s.e.)", m / (sd / std::sqrt(n)), 4);
  }

  // ---- The negative binomial on overdispersed sites.
  {
    calx::random::Stream st(2004, 4, 0, 0);
    const std::size_t S = 1000, T = 5;
    const double r = 2;
    std::vector<double> lam(S, 3.0), p(S * T, 0.5);
    const auto sim = pd::simulate_visits(st, pd::Mixing::negative_binomial, lam, p, T, r);
    const auto I = pd::intercept(S), J = pd::intercept(S * T);
    const auto nb = pd::n_mixture(sim.y, I, J, 120, pd::Mixing::negative_binomial), po = pd::n_mixture(sim.y, I, J, 120);
    s.require("n-mixture (negative binomial): converges, with a positive definite information", nb.converged && nb.information_positive,
              std::to_string(nb.iterations) + " iterations");
    s.check("n-mixture (negative binomial): (log r̂ − log r)/se, on 1,000 sites", (nb.log_size - std::log(r)) / nb.log_size_se, 3);
    s.check("n-mixture (negative binomial): (log λ̂ − log λ)/se", (nb.lambda_coefficients[0] - std::log(3.0)) / nb.lambda_se[0], 3);
    s.require("n-mixture (negative binomial): AIC prefers it to the Poisson", nb.aic < po.aic,
              "ΔAIC " + calx::verification::exponential(po.aic - nb.aic));
  }
  return s.done();
}
