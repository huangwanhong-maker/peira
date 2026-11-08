// Instrumental variables: two-stage least squares, the first stage, and tests of the
// overidentifying restrictions (moved from nomos in peira's stage 2; peira's decision 0004).
//
// THE MODEL. y = W γ + E δ + u with W exogenous (an intercept among its columns, as the caller
// chooses) and E endogenous (correlated with u); the excluded instruments Z_x are correlated with
// E and not with u. Write X = [W E] (K columns, coefficients β = (γ, δ)) and Z = [W Z_x] (L
// columns); identification needs L ≥ K (order condition) and ZᵀX of rank K.
//
// TWO-STAGE LEAST SQUARES (Theil 1953; Basmann 1957): X̂ = P_Z X (W reproduces itself; each column
// of E is replaced by its fitted values on Z, computed by tyche's OLS), then β̂ = (X̂ᵀX̂)⁻¹ X̂ᵀ y =
// (XᵀP_Z X)⁻¹ XᵀP_Z y. When L = K it is the simple instrumental-variables estimator (ZᵀX)⁻¹ Zᵀy
// (checked). The residuals are û = y − Xβ̂ (not y − X̂β̂). Covariances (Wooldridge 2010, ch. 5;
// Cameron and Trivedi 2005, ch. 6), with B = (X̂ᵀX̂)⁻¹:
// - classical: σ̂² B, σ̂² = ûᵀû / (n − K);
// - hc0: B (Σ û_i² x̂_i x̂_iᵀ) B (White 1982); hc1: hc0 · n/(n − K);
// - cluster: G/(G − 1) · (n − 1)/(n − K) · B (Σ_g s_g s_gᵀ) B, s_g = Σ_{i∈g} û_i x̂_i (CR1).
// t statistics on n − K degrees of freedom, or G − 1 with clusters, as tyche's OLS does.
//
// THE FIRST STAGE. For each endogenous regressor, the regression on Z and the F statistic of the
// excluded instruments (the Wald form under the chosen covariance, the classical F when it is
// classical), with its partial R² (Bound, Jaeger and Baker 1995). Weak instruments bias 2SLS
// towards OLS and distort its tests: with one endogenous regressor, Staiger and Stock (1997)
// suggest F > 10, and Stock and Yogo (2005) tabulate the critical values of the (Cragg–Donald)
// statistic for a given maximal bias or size distortion; their tables are not reproduced here.
//
// OVERIDENTIFICATION (L > K). Sargan's (1958) statistic n ûᵀP_Z û / ûᵀû, n times the uncentred R²
// of û on Z, valid under homoskedasticity, and Hansen's (1982) J, robust: with
// Ω̂ = (1/n) Σ û_i² z_i z_iᵀ from the 2SLS residuals (or the cluster sums), the two-step efficient
// GMM estimator β̃ = (XᵀZ Ω̂⁻¹ ZᵀX)⁻¹ XᵀZ Ω̂⁻¹ Zᵀy and J = n ḡ(β̃)ᵀ Ω̂⁻¹ ḡ(β̃), ḡ(β) = Zᵀ(y − Xβ)/n
// (Hayashi 2000, ch. 3; as ivreg2 computes it). Both are χ² with L − K degrees of freedom under
// the null that every instrument is valid; NaN when L = K.

#pragma once

#include <tyche/probability/distributions.hpp>
#include <tyche/regression/linear.hpp>

#include <calx/linalg/cholesky.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <vector>

#include "common.hpp"

namespace peira::quasi {

enum class IvCovariance { classical, hc0, hc1, cluster };

struct IvOptions {
  IvCovariance covariance = IvCovariance::hc1;
  std::vector<std::int64_t> clusters;  // one label per observation, for IvCovariance::cluster
  double confidence = 0.95;
};

struct FirstStage {
  double f_statistic = 0, df1 = 0, df2 = 0, p_value = 0;  // the excluded instruments jointly zero
  double partial_r_squared = 0;
};

struct IvFit {
  std::size_t n = 0, k = 0, instruments = 0;          // n, K = |W| + |E|, L = |W| + |Z_x|
  std::vector<double> beta, se, t, p, ci_lower, ci_upper;  // W's columns, then E's
  Matrix covariance;
  double df_inference = 0;
  std::size_t clusters = 0;
  std::vector<double> residuals;                      // y − Xβ̂
  std::vector<FirstStage> first_stage;                // one per endogenous regressor
  double overidentification_df = 0;                   // L − K
  double sargan = 0, sargan_p_value = 0;
  double hansen_j = 0, hansen_p_value = 0;
  std::vector<double> gmm_beta;                       // the two-step efficient GMM estimate (L > K)
};

inline IvFit two_stage_least_squares(std::span<const double> y, const Rectangular& exogenous, const Rectangular& endogenous, const Rectangular& instruments,
                                     const IvOptions& options = {}) {
  const std::size_t n = y.size();
  const std::size_t kw = exogenous.cols, ke = endogenous.cols, lx = instruments.cols;
  detail::require((kw == 0 || exogenous.rows == n) && (ke == 0 || endogenous.rows == n) && instruments.rows == n, "one row per observation");
  detail::require(lx >= ke, "at least as many excluded instruments as endogenous regressors");
  detail::require(options.covariance != IvCovariance::cluster || options.clusters.size() == n, "a cluster label per observation");
  const Rectangular Z = detail::join_columns(exogenous, instruments), X = detail::join_columns(exogenous, endogenous);
  const std::size_t K = kw + ke, L = kw + lx;
  IvFit fit;
  fit.n = n;
  fit.k = K;
  fit.instruments = L;
  // First stages.
  Rectangular Xh = X;
  tr::OlsOptions fo;
  if (options.covariance == IvCovariance::hc0) fo.covariance = tr::Covariance::hc0;
  if (options.covariance == IvCovariance::hc1) fo.covariance = tr::Covariance::hc1;
  if (options.covariance == IvCovariance::cluster) fo.covariance = tr::Covariance::cluster, fo.clusters = options.clusters;
  for (std::size_t e = 0; e < ke; e++) {
    const auto col = endogenous.column(e);
    const auto first = tr::ols(Z, col, fo);
    for (std::size_t i = 0; i < n; i++) Xh(i, kw + e) = first.fitted[i];
    Rectangular R(lx, L);
    for (std::size_t a = 0; a < lx; a++) R(a, kw + a) = 1;
    const auto w = tr::wald(first, R, std::vector<double>(lx, 0.0));
    FirstStage fs;
    fs.f_statistic = w.f_statistic;
    fs.df1 = w.df1;
    fs.df2 = w.df2;
    fs.p_value = w.f_p_value;
    double restricted = 0;
    if (kw > 0) {
      tr::OlsOptions plain;
      plain.drop_collinear = true;
      restricted = tr::ols(exogenous, col, plain).rss;
    } else {
      for (double v : col) restricted += v * v;
    }
    fs.partial_r_squared = 1 - first.rss / restricted;
    fit.first_stage.push_back(fs);
  }
  // The second stage.
  const auto second = tr::ols(Xh, y);
  fit.beta = second.beta;
  fit.residuals.resize(n);
  double uu = 0;
  for (std::size_t i = 0; i < n; i++) {
    double xb = 0;
    for (std::size_t j = 0; j < K; j++) xb += X(i, j) * fit.beta[j];
    fit.residuals[i] = y[i] - xb;
    uu += fit.residuals[i] * fit.residuals[i];
  }
  const Matrix& B = second.xtx_inverse;
  const double dn = static_cast<double>(n), dk = static_cast<double>(K);
  std::size_t groups = 0;
  switch (options.covariance) {
    case IvCovariance::classical: fit.covariance = (uu / (dn - dk)) * B; break;
    case IvCovariance::hc0: fit.covariance = detail::sandwich(B, detail::sandwich_meat(Xh, fit.residuals, {}, groups)); break;
    case IvCovariance::hc1: fit.covariance = (dn / (dn - dk)) * detail::sandwich(B, detail::sandwich_meat(Xh, fit.residuals, {}, groups)); break;
    case IvCovariance::cluster: {
      const Matrix M = detail::sandwich_meat(Xh, fit.residuals, options.clusters, groups);
      detail::require(groups >= 2, "at least two clusters");
      const double G = static_cast<double>(groups);
      fit.covariance = (G / (G - 1) * (dn - 1) / (dn - dk)) * detail::sandwich(B, M);
      fit.clusters = groups;
      break;
    }
  }
  fit.df_inference = options.covariance == IvCovariance::cluster ? static_cast<double>(groups) - 1 : dn - dk;
  const tyche::probability::StudentT tdist{fit.df_inference};
  const double tcrit = tdist.upper_quantile((1 - options.confidence) / 2);
  for (std::size_t j = 0; j < K; j++) {
    const double se = std::sqrt(fit.covariance(j, j));
    fit.se.push_back(se);
    fit.t.push_back(fit.beta[j] / se);
    fit.p.push_back(std::min(1.0, 2 * tdist.survival(std::fabs(fit.beta[j] / se))));
    fit.ci_lower.push_back(fit.beta[j] - tcrit * se);
    fit.ci_upper.push_back(fit.beta[j] + tcrit * se);
  }
  // Overidentification.
  fit.overidentification_df = static_cast<double>(L - K);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (L == K) {
    fit.sargan = fit.sargan_p_value = fit.hansen_j = fit.hansen_p_value = nan;
    return fit;
  }
  const auto aux = tr::ols(Z, fit.residuals);
  double explained = 0;
  for (double v : aux.fitted) explained += v * v;
  fit.sargan = dn * explained / uu;
  const tyche::probability::ChiSquared chi{fit.overidentification_df};
  fit.sargan_p_value = chi.survival(fit.sargan);
  // Hansen's J at the two-step efficient GMM estimate.
  std::size_t g2 = 0;
  Matrix Omega = (1 / dn) * detail::sandwich_meat(Z, fit.residuals, options.covariance == IvCovariance::cluster ? std::span<const std::int64_t>(options.clusters) : std::span<const std::int64_t>{}, g2);
  // ZᵀX (L × K) and Zᵀy.
  Rectangular ZX(L, K);
  std::vector<double> Zy(L, 0.0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t a = 0; a < L; a++) {
      Zy[a] += Z(i, a) * y[i];
      for (std::size_t b = 0; b < K; b++) ZX(a, b) += Z(i, a) * X(i, b);
    }
  // Ω⁻¹ ZᵀX column by column, and Ω⁻¹ Zᵀy.
  Rectangular OZX(L, K);
  for (std::size_t b = 0; b < K; b++) {
    const auto c = detail::solve_spd(Omega, ZX.column(b));
    for (std::size_t a = 0; a < L; a++) OZX(a, b) = c[a];
  }
  const auto OZy = detail::solve_spd(Omega, Zy);
  Matrix A(K);
  std::vector<double> rhs(K, 0.0);
  for (std::size_t a = 0; a < K; a++) {
    for (std::size_t b = 0; b < K; b++)
      for (std::size_t l = 0; l < L; l++) A(a, b) += ZX(l, a) * OZX(l, b);
    for (std::size_t l = 0; l < L; l++) rhs[a] += ZX(l, a) * OZy[l];
  }
  fit.gmm_beta = detail::solve_spd(A, rhs);
  std::vector<double> g(L, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    double r = y[i];
    for (std::size_t b = 0; b < K; b++) r -= X(i, b) * fit.gmm_beta[b];
    for (std::size_t a = 0; a < L; a++) g[a] += Z(i, a) * r / dn;
  }
  const auto Og = detail::solve_spd(Omega, g);
  fit.hansen_j = 0;
  for (std::size_t a = 0; a < L; a++) fit.hansen_j += dn * g[a] * Og[a];
  fit.hansen_p_value = chi.survival(fit.hansen_j);
  return fit;
}

// THE WALD ESTIMATOR. With a binary instrument Z (an offer, an encouragement, a lottery) and a binary treatment d, under
// independence, exclusion and monotonicity (no defiers), the ratio
//   τ̂ = (ȳ₁ − ȳ₀) / (d̄₁ − d̄₀),
// means among Z = 1 and Z = 0, estimates the local average treatment effect: the average effect on compliers, those whom the
// instrument moves (Wald 1940; Imbens and Angrist 1994; Angrist, Imbens and Rubin 1996). It is 2SLS of y on an intercept and d
// with Z the instrument, exactly (checked), whose robust covariance gives its standard error; the first stage's difference d̄₁ − d̄₀
// is the compliers' share.
struct Wald {
  double estimate = 0, se = 0, t = 0, p = 0, ci_lower = 0, ci_upper = 0;
  double compliers = 0;                 // d̄₁ − d̄₀
  double reduced_form = 0;              // ȳ₁ − ȳ₀
  IvFit fit;
};
inline Wald wald(std::span<const double> y, std::span<const double> d, std::span<const int> z, const IvOptions& options = {}) {
  const std::size_t n = y.size();
  detail::require(d.size() == n && z.size() == n, "a treatment and an instrument for every observation");
  double y1 = 0, y0 = 0, d1 = 0, d0 = 0, n1 = 0, n0 = 0;
  Rectangular W(n, 1), E(n, 1), Z(n, 1);
  for (std::size_t i = 0; i < n; i++) {
    detail::require(z[i] == 0 || z[i] == 1, "an instrument of 0 or 1");
    (z[i] ? y1 : y0) += y[i], (z[i] ? d1 : d0) += d[i], (z[i] ? n1 : n0) += 1;
    W(i, 0) = 1, E(i, 0) = d[i], Z(i, 0) = z[i];
  }
  detail::require(n1 > 0 && n0 > 0, "observations with each value of the instrument");
  Wald w;
  w.reduced_form = y1 / n1 - y0 / n0;
  w.compliers = d1 / n1 - d0 / n0;
  detail::require(std::fabs(w.compliers) > 0, "an instrument that moves the treatment");
  w.fit = two_stage_least_squares(y, W, E, Z, options);
  w.estimate = w.fit.beta[1], w.se = w.fit.se[1], w.t = w.fit.t[1], w.p = w.fit.p[1], w.ci_lower = w.fit.ci_lower[1], w.ci_upper = w.fit.ci_upper[1];
  return w;
}

}  // namespace peira::quasi
