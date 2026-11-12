// Latent counts seen through binomial detection: the likelihood that occupancy and the N-mixture model share (docs/theory.md,
// §10.2 and §10.3).
//
// THE MODEL. At site i a latent count N_i has a law whose mean depends on x_i'α (x_i the site's row of the design X):
// - bernoulli: N_i ∈ {0, 1}, P(N_i = 1) = ψ_i = logistic(x_i'α), the site's occupancy (MacKenzie et al. 2002);
// - poisson: N_i ~ Poisson(λ_i), λ_i = exp(x_i'α), its abundance (Royle 2004);
// - negative_binomial: N_i with mean λ_i = exp(x_i'α), size r = exp(ρ) and variance λ_i + λ_i²/r (Royle 2004), which tends
//   to the Poisson as r → ∞.
// Given N_i, the count on each visit is y_ij ~ Binomial(N_i, p_ij), p_ij = logistic(w_ij'β), independently over visits (w_ij
// the row i·visits + j of the design W). The site's likelihood sums over what it holds:
//   L_i = Σ_{N = max_j y_ij}^{K} P(N_i = N) Π_j C(N, y_ij) p_ij^{y_ij} (1 − p_ij)^{N − y_ij},
// with K = 1 for presence and, for abundance, a bound the caller states: the true likelihood is the limit K → ∞, and the bound's
// effect is reported by the N-mixture fit. Visits not made (y = −1) are left out of the product; a site never visited carries
// no information and is left out of the likelihood.
//
// DERIVATIVES. The gradient is the posterior expectation of the complete-data score, and the Hessian is Louis's identity
// (Louis 1982, "Finding the observed information matrix when using the EM algorithm", J. R. Statist. Soc. B 44, 226–233):
//   ∇²ℓ_i = E[∇² log f(y_i, N) | y_i] + Var[∇ log f(y_i, N) | y_i],
// computed exactly by the same sum over N. The complete-data log density separates into the latent law's part (in x_i'α and ρ)
// and detection's (in each w_ij'β), and detection's score in w_ij'β is y_ij − N p_ij, linear in N; so the sums need only the
// posterior moments of N, of the latent law's scores and of their products, taken about their means. Each site's sum runs in
// log space (log-sum-exp), so that counts of hundreds lose nothing to underflow.
//
// For the negative binomial, Γ(N + r)/Γ(r) = Π_{k<N} (r + k), ψ(N + r) − ψ(r) = Σ_{k<N} 1/(r + k) and
// ψ'(N + r) − ψ'(r) = −Σ_{k<N} 1/(r + k)², tabulated once per evaluation, exact to rounding for every N ≤ K.
//
// FITTING. calx's BFGS on the mean negative log-likelihood per site with this gradient, from the caller's start (one restart
// from where it stopped if it stopped unconverged), then Newton's steps with the exact Hessian while they reduce the gradient,
// so that the estimate is found to rounding (as survival's Cox fit does). Converged means a mean gradient per site below the
// tolerance. The covariance is the inverse of the observed information at the estimate.

#pragma once

#include <calx/optimisation/quasi_newton.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include "common.hpp"

namespace peira::detection {

enum class Mixing { bernoulli, poisson, negative_binomial };

namespace detail {

inline std::size_t latent_parameters(const Rectangular& X, const Rectangular& W, Mixing m) {
  return X.cols + W.cols + (m == Mixing::negative_binomial ? 1 : 0);
}

inline void latent_validate(const SiteVisits& y, const Rectangular& X, const Rectangular& W, Mixing m, int K) {
  require(y.y.size() == y.sites * y.visits, "a table of sites by visits");
  require(X.rows == y.sites, "one row of the latent law's design per site");
  require(W.rows == y.sites * y.visits, "one row of detection's design per site and visit (row i·visits + j)");
  require(K >= 1, "an upper bound K ≥ 1 for the latent count");
  for (int v : y.y) {
    require(v >= -1, "counts ≥ 0, or −1 for a visit not made");
    require(v <= K, "every count at most the bound K");
    if (m == Mixing::bernoulli) require(v <= 1, "detections 0 or 1 for presence");
  }
}

struct Latent {
  double loglik = 0;
  std::vector<double> gradient;
  Matrix hessian;
  std::size_t sites_used = 0;
};

// ℓ(θ), and with order ≥ 1 its gradient, with order 2 its Hessian; θ = (α, β[, ρ]). With `posterior`, each site's law of N given
// its counts (the prior for a site never visited), N = 0 … K.
inline Latent latent_evaluate(const SiteVisits& y, const Rectangular& X, const Rectangular& W, Mixing mixing, int K, std::span<const double> theta,
                              int order, std::vector<std::vector<double>>* posterior = nullptr) {
  const std::size_t kx = X.cols, kw = W.cols, k = latent_parameters(X, W, mixing), T = y.visits;
  require(theta.size() == k, "one coefficient per column of the designs (and log r for the negative binomial)");
  const bool nb = mixing == Mixing::negative_binomial;
  const std::size_t KK = static_cast<std::size_t>(K);
  Latent out;
  out.gradient.assign(k, 0.0);
  out.hessian = Matrix(k);
  if (posterior) posterior->assign(y.sites, std::vector<double>(KK + 1, 0.0));

  std::vector<double> lf(KK + 1);
  for (std::size_t n = 0; n <= KK; n++) lf[n] = calx::special_functions::log_gamma(static_cast<double>(n) + 1);
  const double rho = nb ? theta[k - 1] : 0, r = nb ? std::exp(rho) : 0;
  std::vector<double> G, D, Tr;  // Σ_{k<N} log(r + k), 1/(r + k), 1/(r + k)²
  if (nb) {
    G.assign(KK + 1, 0.0), D.assign(KK + 1, 0.0), Tr.assign(KK + 1, 0.0);
    for (std::size_t n = 1; n <= KK; n++) {
      const double a = r + static_cast<double>(n - 1);
      G[n] = G[n - 1] + std::log(a), D[n] = D[n - 1] + 1 / a, Tr[n] = Tr[n - 1] + 1 / (a * a);
    }
  }

  std::vector<int> yy;
  std::vector<std::size_t> jj;
  std::vector<double> p, lj(KK + 1), wN(KK + 1), ss(KK + 1), sr(KK + 1), hs(KK + 1), hsr(KK + 1), hrr(KK + 1);
  for (std::size_t i = 0; i < y.sites; i++) {
    double eta = 0;
    for (std::size_t a = 0; a < kx; a++) eta += X(i, a) * theta[a];
    // The latent law at N = 0 … K: log P(N), its scores in η = x'α (and ρ), and their second derivatives.
    const double psi = logistic(eta), lambda = std::exp(eta);
    auto log_prior = [&](std::size_t n) {
      const double nd = static_cast<double>(n);
      switch (mixing) {
        case Mixing::bernoulli: return n == 0 ? -log1pexp(eta) : -log1pexp(-eta);
        case Mixing::poisson: return -lambda + nd * eta - lf[n];
        default: return G[n] - lf[n] - r * std::log1p(lambda / r) + nd * (eta - std::log(r + lambda));
      }
    };
    yy.clear(), jj.clear(), p.clear();
    int ymax = 0;
    double A = 0, B = 0, C = 0;  // Σ y η_d, Σ log(1 − p), −Σ log y!
    for (std::size_t j = 0; j < T; j++) {
      const int v = y(i, j);
      if (v < 0) continue;
      double ed = 0;
      for (std::size_t b = 0; b < kw; b++) ed += W(i * T + j, b) * theta[kx + b];
      yy.push_back(v), jj.push_back(j), p.push_back(logistic(ed));
      A += v * ed, B -= log1pexp(ed), C -= lf[static_cast<std::size_t>(v)];
      ymax = std::max(ymax, v);
    }
    if (yy.empty()) {
      if (posterior) {
        double s = 0;
        for (std::size_t n = 0; n <= (mixing == Mixing::bernoulli ? 1 : KK); n++) s += ((*posterior)[i][n] = std::exp(log_prior(n)));
        for (double& v : (*posterior)[i]) v /= s;
      }
      continue;
    }
    out.sites_used++;
    const std::size_t lo = static_cast<std::size_t>(ymax), hi = mixing == Mixing::bernoulli ? 1 : KK;
    const double J = static_cast<double>(yy.size());
    double m = -std::numeric_limits<double>::infinity();
    for (std::size_t n = lo; n <= hi; n++) {
      double l = log_prior(n) + A + static_cast<double>(n) * B + C + J * lf[n];
      for (int v : yy) l -= lf[n - static_cast<std::size_t>(v)];
      lj[n] = l;
      m = std::max(m, l);
    }
    double s = 0;
    for (std::size_t n = lo; n <= hi; n++) s += (wN[n] = std::exp(lj[n] - m));
    out.loglik += m + std::log(s);
    for (std::size_t n = lo; n <= hi; n++) wN[n] /= s;
    if (posterior)
      for (std::size_t n = lo; n <= hi; n++) (*posterior)[i][n] = wN[n];
    if (order < 1) continue;

    // The scores and second derivatives of log P(N) in η (and ρ).
    double EN = 0, Es = 0, Er = 0, Eh = 0, Ehsr = 0, Ehrr = 0;
    for (std::size_t n = lo; n <= hi; n++) {
      const double nd = static_cast<double>(n);
      switch (mixing) {
        case Mixing::bernoulli: ss[n] = nd - psi, hs[n] = -psi * (1 - psi); break;
        case Mixing::poisson: ss[n] = nd - lambda, hs[n] = -lambda; break;
        default: {
          const double rl = r + lambda, g = D[n] - std::log1p(lambda / r) + (lambda - nd) / rl;
          const double gp = -Tr[n] + 1 / r - 1 / rl - (lambda - nd) / (rl * rl);
          ss[n] = r * (nd - lambda) / rl, hs[n] = -r * lambda * (r + nd) / (rl * rl);
          sr[n] = r * g, hrr[n] = r * g + r * r * gp, hsr[n] = r * lambda * (nd - lambda) / (rl * rl);
        }
      }
      EN += wN[n] * nd, Es += wN[n] * ss[n], Eh += wN[n] * hs[n];
      if (nb) Er += wN[n] * sr[n], Ehsr += wN[n] * hsr[n], Ehrr += wN[n] * hrr[n];
    }
    for (std::size_t a = 0; a < kx; a++) out.gradient[a] += Es * X(i, a);
    for (std::size_t q = 0; q < yy.size(); q++)
      for (std::size_t b = 0; b < kw; b++) out.gradient[kx + b] += (yy[q] - EN * p[q]) * W(i * T + jj[q], b);
    if (nb) out.gradient[k - 1] += Er;
    if (order < 2) continue;

    double VN = 0, Vs = 0, CsN = 0, Vr = 0, Csr = 0, CrN = 0;  // posterior (co)variances, about the means
    for (std::size_t n = lo; n <= hi; n++) {
      const double dN = static_cast<double>(n) - EN, ds = ss[n] - Es, dr = nb ? sr[n] - Er : 0;
      VN += wN[n] * dN * dN, Vs += wN[n] * ds * ds, CsN += wN[n] * ds * dN;
      if (nb) Vr += wN[n] * dr * dr, Csr += wN[n] * ds * dr, CrN += wN[n] * dr * dN;
    }
    std::vector<double> bvec(kw, 0.0);  // Σ_j p_ij w_ij
    for (std::size_t q = 0; q < yy.size(); q++) {
      const std::size_t row = i * T + jj[q];
      for (std::size_t b = 0; b < kw; b++) bvec[b] += p[q] * W(row, b);
      for (std::size_t b = 0; b < kw; b++)
        for (std::size_t c = 0; c < kw; c++) out.hessian(kx + b, kx + c) -= EN * p[q] * (1 - p[q]) * W(row, b) * W(row, c);
    }
    for (std::size_t a = 0; a < kx; a++) {
      for (std::size_t b = 0; b < kx; b++) out.hessian(a, b) += (Eh + Vs) * X(i, a) * X(i, b);
      for (std::size_t b = 0; b < kw; b++) {
        out.hessian(a, kx + b) -= CsN * X(i, a) * bvec[b];
        out.hessian(kx + b, a) -= CsN * X(i, a) * bvec[b];
      }
      if (nb) {
        out.hessian(a, k - 1) += (Ehsr + Csr) * X(i, a);
        out.hessian(k - 1, a) += (Ehsr + Csr) * X(i, a);
      }
    }
    for (std::size_t b = 0; b < kw; b++) {
      for (std::size_t c = 0; c < kw; c++) out.hessian(kx + b, kx + c) += VN * bvec[b] * bvec[c];
      if (nb) {
        out.hessian(kx + b, k - 1) -= CrN * bvec[b];
        out.hessian(k - 1, kx + b) -= CrN * bvec[b];
      }
    }
    if (nb) out.hessian(k - 1, k - 1) += Ehrr + Vr;
  }
  return out;
}

struct LatentFit {
  std::vector<double> theta, se;
  Matrix covariance;
  double loglik = 0, gradient_norm = 0;
  bool converged = false, information_positive = false;
  int iterations = 0;
  std::size_t sites_used = 0;
};

inline LatentFit latent_fit(const SiteVisits& y, const Rectangular& X, const Rectangular& W, Mixing mixing, int K, std::vector<double> theta0,
                            const FitOptions& o) {
  const std::size_t k = latent_parameters(X, W, mixing);
  require(theta0.size() == k, "a start with one coefficient per parameter");
  const auto first = latent_evaluate(y, X, W, mixing, K, theta0, 0);
  require(first.sites_used > 0, "at least one site visited");
  const double scale = 1.0 / static_cast<double>(first.sites_used);
  auto f = [&](std::span<const double> t) {
    const double l = latent_evaluate(y, X, W, mixing, K, t, 0).loglik;
    return std::isfinite(l) ? -l * scale : std::numeric_limits<double>::infinity();
  };
  auto g = [&](std::span<const double> t) {
    auto e = latent_evaluate(y, X, W, mixing, K, t, 1);
    for (double& v : e.gradient) v *= -scale;
    return e.gradient;
  };
  calx::optimisation::BfgsOptions bo;
  bo.gradient_tolerance = o.gradient_tolerance;
  bo.max_iterations = o.max_iterations;
  auto res = calx::optimisation::bfgs(f, g, std::move(theta0), bo);
  int iterations = res.iterations;
  if (!res.converged) {
    res = calx::optimisation::bfgs(f, g, res.x, bo);
    iterations += res.iterations;
  }
  LatentFit fit;
  fit.theta = res.x;
  fit.iterations = iterations;
  // Newton's steps with the exact Hessian, from where BFGS stopped: its line search compares values, which near the maximum
  // differ by less than their rounding, while Newton's step is computed from the gradient and keeps converging. A step is
  // taken only where the information is positive definite and the gradient falls.
  auto e = latent_evaluate(y, X, W, mixing, K, fit.theta, 2);
  auto norm = [&](const std::vector<double>& v) {
    double m = 0;
    for (double a : v) m = std::max(m, std::fabs(a) * scale);
    return m;
  };
  for (int step = 0; step < 20 && norm(e.gradient) > 1e-14; step++) {
    Matrix info(k);
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t b = 0; b < k; b++) info(a, b) = -e.hessian(a, b);
    const auto c = calx::linalg::cholesky(info);
    if (!c.positive_definite) break;
    const auto dx = calx::linalg::cholesky_solve(c, e.gradient);
    std::vector<double> next = fit.theta;
    for (std::size_t a = 0; a < k; a++) next[a] += dx[a];
    auto en = latent_evaluate(y, X, W, mixing, K, next, 2);
    if (!std::isfinite(en.loglik) || norm(en.gradient) >= norm(e.gradient)) break;
    fit.theta = std::move(next), e = std::move(en);
    fit.iterations++;
  }
  fit.gradient_norm = norm(e.gradient);
  fit.converged = fit.gradient_norm <= std::max(o.gradient_tolerance, 1e-12) || res.converged;
  fit.loglik = e.loglik;
  fit.sites_used = e.sites_used;
  fit.information_positive = information_inverse(e.hessian, fit.covariance);
  if (!fit.information_positive) fit.covariance = Matrix(k);
  fit.se = standard_errors(fit.covariance, fit.information_positive, k);
  return fit;
}

// A start for the latent law's coefficients: on each column of X that is one constant c ≠ 0 at every site (an intercept), the
// value v/c, so that the linear predictor starts at v; elsewhere 0.
inline void start_on_constant_columns(const Rectangular& X, double v, std::vector<double>& theta, std::size_t offset) {
  for (std::size_t a = 0; a < X.cols; a++) {
    if (X.rows == 0) break;
    const double c = X(0, a);
    bool constant = c != 0;
    for (std::size_t i = 1; i < X.rows && constant; i++) constant = X(i, a) == c;
    if (constant) {
      theta[offset + a] = v / c;
      break;  // one intercept is enough
    }
  }
}

}  // namespace detail

}  // namespace peira::detection
