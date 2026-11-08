// Times to an event, seen through censoring (docs/theory.md, §5).
//
// Each unit has a time T_i to its event and a time C_i at which it leaves observation; what is seen is X_i = min(T_i, C_i)
// and δ_i = [T_i ≤ C_i]. Censoring is a coarsening of the full data: under independent censoring the hazard of the
// event among those still at risk is identified at every time, and everything here is built from the risk sets
// R(t) = {i : X_i ≥ t} and the events at t.
//
// - kaplan_meier: S(t) = Π_{t_j ≤ t} (1 − d_j/n_j) (Kaplan and Meier 1958), Greenwood's (1926) variance
//   S² Σ d_j/(n_j(n_j − d_j)), and the log(−log) interval S^{exp(±z σ)}, σ² = Σ d/(n(n − d)) / (log S)²
//   (Kalbfleisch and Prentice 2002, The Statistical Analysis of Failure Time Data, 2nd ed., §1.4); the median, the first
//   time S ≤ 1/2.
// - nelson_aalen: H(t) = Σ_{t_j ≤ t} d_j/n_j (Nelson 1972, Aalen 1978), with variance Σ d_j/n_j².
// - log_rank: k groups compared by O − E at each event time, with the hypergeometric covariance; χ² on k − 1 degrees of
//   freedom (Mantel 1966; Peto and Peto 1972).
// - cox: the proportional hazards model λ(t | x) = λ₀(t) exp(xβ) (Cox 1972), β̂ maximising the partial likelihood with
//   Breslow's (1974) or Efron's (1977) handling of tied times, by Newton–Raphson with step halving; standard errors from
//   the inverse of the observed information; Breslow's estimate of the baseline cumulative hazard.
// - aalen_johansen: with competing causes, the cumulative incidence of cause c, F_c(t) = Σ_{t_j ≤ t} S(t_j−) d_cj/n_j,
//   S the Kaplan–Meier estimate of surviving every cause (Aalen and Johansen 1978); the incidences sum to 1 − S exactly.
//
// Limits: censoring is assumed independent of the event given what the model conditions on; Cox's model assumes its
// hazards proportional, which is not checked here.

#pragma once

#include <tyche/probability/distributions.hpp>

#include <calx/linalg/cholesky.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/linalg/rectangular.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace peira::survival {

using calx::linalg::Matrix;
using calx::linalg::Rectangular;

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::survival: ") + what);
}
inline double z_two_sided(double confidence) { return tyche::probability::Normal{0, 1}.upper_quantile((1 - confidence) / 2); }

// The distinct times at which events occur, with the numbers at risk and the events (in total, and by a label).
struct Table {
  std::vector<double> time, at_risk, events;
};
inline Table table(std::span<const double> time, std::span<const int> event) {
  require(time.size() == event.size(), "one event indicator per time");
  std::map<double, std::pair<double, double>> by_time;  // time → (leaving at it, events at it)
  for (std::size_t i = 0; i < time.size(); i++) {
    auto& e = by_time[time[i]];
    e.first++;
    e.second += event[i] ? 1 : 0;
  }
  Table t;
  double at_risk = static_cast<double>(time.size());
  for (const auto& [x, e] : by_time) {
    if (e.second > 0) { t.time.push_back(x); t.at_risk.push_back(at_risk); t.events.push_back(e.second); }
    at_risk -= e.first;
  }
  return t;
}

inline Matrix inverse(const Matrix& A) {
  const auto f = calx::linalg::cholesky(A);
  require(f.positive_definite, "a positive definite information matrix");
  Matrix inv(A.n);
  std::vector<double> e(A.n, 0.0);
  for (std::size_t j = 0; j < A.n; j++) {
    std::fill(e.begin(), e.end(), 0.0);
    e[j] = 1;
    const auto col = calx::linalg::cholesky_solve(f, e);
    for (std::size_t i = 0; i < A.n; i++) inv(i, j) = col[i];
  }
  return inv;
}
}  // namespace detail

struct KaplanMeier {
  std::vector<double> time, at_risk, events, survival, se, lower, upper;
  double median = std::numeric_limits<double>::quiet_NaN();  // NaN when S never falls to 1/2
};

inline KaplanMeier kaplan_meier(std::span<const double> time, std::span<const int> event, double confidence = 0.95) {
  const auto t = detail::table(time, event);
  KaplanMeier km;
  km.time = t.time, km.at_risk = t.at_risk, km.events = t.events;
  const double z = detail::z_two_sided(confidence);
  double S = 1, g = 0;  // g = Σ d/(n(n − d))
  for (std::size_t j = 0; j < t.time.size(); j++) {
    const double n = t.at_risk[j], d = t.events[j];
    S *= 1 - d / n;
    g += n > d ? d / (n * (n - d)) : std::numeric_limits<double>::infinity();
    km.survival.push_back(S);
    km.se.push_back(S * std::sqrt(g));
    if (S > 0 && S < 1 && std::isfinite(g)) {
      const double sigma = std::sqrt(g) / std::fabs(std::log(S));
      km.lower.push_back(std::pow(S, std::exp(z * sigma)));
      km.upper.push_back(std::pow(S, std::exp(-z * sigma)));
    } else {
      km.lower.push_back(std::numeric_limits<double>::quiet_NaN());
      km.upper.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    if (std::isnan(km.median) && S <= 0.5) km.median = t.time[j];
  }
  return km;
}

// S at time x, from an estimate (a step function, right-continuous; 1 before the first event).
inline double survival_at(const KaplanMeier& km, double x) {
  double S = 1;
  for (std::size_t j = 0; j < km.time.size() && km.time[j] <= x; j++) S = km.survival[j];
  return S;
}

struct NelsonAalen {
  std::vector<double> time, hazard, variance;
};

inline NelsonAalen nelson_aalen(std::span<const double> time, std::span<const int> event) {
  const auto t = detail::table(time, event);
  NelsonAalen na;
  na.time = t.time;
  double H = 0, v = 0;
  for (std::size_t j = 0; j < t.time.size(); j++) {
    H += t.events[j] / t.at_risk[j];
    v += t.events[j] / (t.at_risk[j] * t.at_risk[j]);
    na.hazard.push_back(H), na.variance.push_back(v);
  }
  return na;
}

struct LogRank {
  double statistic = 0, df = 0, p_value = 1;
  std::vector<int> groups;                 // the labels, in increasing order
  std::vector<double> observed, expected;  // per group
};

inline LogRank log_rank(std::span<const double> time, std::span<const int> event, std::span<const int> group) {
  detail::require(time.size() == event.size() && time.size() == group.size(), "one event and one group per time");
  std::map<int, std::size_t> label;
  for (int g : group) label.emplace(g, 0);
  std::size_t k = 0;
  for (auto& [g, index] : label) index = k++;
  detail::require(k >= 2, "at least two groups");
  // By time: those leaving, and the events, per group.
  std::map<double, std::pair<std::vector<double>, std::vector<double>>> by_time;
  for (std::size_t i = 0; i < time.size(); i++) {
    auto& e = by_time[time[i]];
    if (e.first.empty()) e.first.assign(k, 0.0), e.second.assign(k, 0.0);
    e.first[label[group[i]]]++;
    if (event[i]) e.second[label[group[i]]]++;
  }
  std::vector<double> at_risk(k, 0.0);
  for (int g : group) at_risk[label[g]]++;
  LogRank r;
  for (const auto& [g, index] : label) r.groups.push_back(g);
  r.observed.assign(k, 0.0), r.expected.assign(k, 0.0);
  Matrix V(k);
  for (const auto& [x, e] : by_time) {
    const double n = std::accumulate(at_risk.begin(), at_risk.end(), 0.0), d = std::accumulate(e.second.begin(), e.second.end(), 0.0);
    if (d > 0) {
      for (std::size_t g = 0; g < k; g++) {
        r.observed[g] += e.second[g];
        r.expected[g] += d * at_risk[g] / n;
      }
      if (n > 1)
        for (std::size_t g = 0; g < k; g++)
          for (std::size_t h = 0; h < k; h++)
            V(g, h) += d * (n - d) / (n - 1) * at_risk[g] / n * ((g == h ? 1.0 : 0.0) - at_risk[h] / n);
    }
    for (std::size_t g = 0; g < k; g++) at_risk[g] -= e.first[g];
  }
  // (O − E)ᵀ V⁻ (O − E) on the first k − 1 groups.
  Matrix W(k - 1);
  std::vector<double> u(k - 1);
  for (std::size_t g = 0; g + 1 < k; g++) {
    u[g] = r.observed[g] - r.expected[g];
    for (std::size_t h = 0; h + 1 < k; h++) W(g, h) = V(g, h);
  }
  const auto f = calx::linalg::cholesky(W);
  detail::require(f.positive_definite, "a non-degenerate covariance of O − E");
  const auto x = calx::linalg::cholesky_solve(f, u);
  for (std::size_t g = 0; g + 1 < k; g++) r.statistic += u[g] * x[g];
  r.df = static_cast<double>(k - 1);
  r.p_value = tyche::probability::ChiSquared{r.df}.survival(r.statistic);
  return r;
}

enum class Ties { breslow, efron };

struct CoxOptions {
  Ties ties = Ties::efron;
  int max_iterations = 50;
  double tolerance = 1e-12;  // on |Δℓ| / (|ℓ| + 1)
  double confidence = 0.95;
};

struct Cox {
  std::vector<double> beta, se, z, p, ci_lower, ci_upper;
  Matrix covariance;
  double log_partial_likelihood = 0, null_log_partial_likelihood = 0;
  int iterations = 0;
  bool converged = false;
  std::vector<double> baseline_time, baseline_hazard;  // Breslow's cumulative baseline hazard at the event times
};

namespace detail {
struct CoxState {
  double loglik = 0;
  std::vector<double> gradient;
  Matrix information;  // −∂²ℓ
};

// ℓ, its gradient and information at β: the distinct times in decreasing order, the risk sums accumulated.
inline CoxState cox_state(const Rectangular& X, std::span<const double> time, std::span<const int> event, std::span<const double> beta, Ties ties) {
  const std::size_t n = X.rows, k = X.cols;
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return time[a] > time[b]; });
  std::vector<double> eta(n);
  for (std::size_t i = 0; i < n; i++) {
    eta[i] = 0;
    for (std::size_t j = 0; j < k; j++) eta[i] += X(i, j) * beta[j];
  }
  CoxState st;
  st.gradient.assign(k, 0.0);
  st.information = Matrix(k);
  double s0 = 0;
  std::vector<double> s1(k, 0.0), s2(k * k, 0.0);
  std::size_t at = 0;
  while (at < n) {
    const double t = time[order[at]];
    // The units at this time join the risk set; its events are gathered apart.
    double d0 = 0, d = 0;
    std::vector<double> d1(k, 0.0), d2(k * k, 0.0), xsum(k, 0.0);
    double eta_sum = 0;
    std::size_t end = at;
    while (end < n && time[order[end]] == t) {
      const std::size_t i = order[end];
      const double w = std::exp(eta[i]);
      s0 += w;
      for (std::size_t a = 0; a < k; a++) {
        s1[a] += w * X(i, a);
        for (std::size_t b = 0; b < k; b++) s2[a * k + b] += w * X(i, a) * X(i, b);
      }
      if (event[i]) {
        d++;
        d0 += w;
        eta_sum += eta[i];
        for (std::size_t a = 0; a < k; a++) {
          d1[a] += w * X(i, a);
          xsum[a] += X(i, a);
          for (std::size_t b = 0; b < k; b++) d2[a * k + b] += w * X(i, a) * X(i, b);
        }
      }
      end++;
    }
    if (d > 0) {
      st.loglik += eta_sum;
      for (std::size_t a = 0; a < k; a++) st.gradient[a] += xsum[a];
      for (int l = 0; l < static_cast<int>(d); l++) {
        const double frac = ties == Ties::efron ? l / d : 0.0;
        const double A = s0 - frac * d0;
        st.loglik -= std::log(A);
        for (std::size_t a = 0; a < k; a++) {
          const double Ba = s1[a] - frac * d1[a];
          st.gradient[a] -= Ba / A;
          for (std::size_t b = 0; b < k; b++) {
            const double Bb = s1[b] - frac * d1[b], C = s2[a * k + b] - frac * d2[a * k + b];
            st.information(a, b) += C / A - Ba * Bb / (A * A);
          }
        }
      }
    }
    at = end;
  }
  return st;
}
}  // namespace detail

// X: one row of covariates per unit, without an intercept (the baseline hazard absorbs it).
inline Cox cox(const Rectangular& X, std::span<const double> time, std::span<const int> event, const CoxOptions& o = {}) {
  const std::size_t n = X.rows, k = X.cols;
  detail::require(time.size() == n && event.size() == n && k >= 1, "one time and one event per row, at least one covariate");
  Cox r;
  std::vector<double> beta(k, 0.0);
  auto st = detail::cox_state(X, time, event, beta, o.ties);
  r.null_log_partial_likelihood = st.loglik;
  for (r.iterations = 1; r.iterations <= o.max_iterations; r.iterations++) {
    const auto f = calx::linalg::cholesky(st.information);
    detail::require(f.positive_definite, "a positive definite information matrix (covariates collinear, or separation)");
    const auto step = calx::linalg::cholesky_solve(f, st.gradient);
    double scale = 1;
    detail::CoxState next;
    std::vector<double> trial(k);
    for (int halving = 0; halving < 30; halving++) {
      for (std::size_t j = 0; j < k; j++) trial[j] = beta[j] + scale * step[j];
      next = detail::cox_state(X, time, event, trial, o.ties);
      if (next.loglik >= st.loglik - 1e-12 * (std::fabs(st.loglik) + 1)) break;
      scale /= 2;
    }
    const double change = std::fabs(next.loglik - st.loglik) / (std::fabs(st.loglik) + 1);
    beta = trial, st = std::move(next);
    if (change < o.tolerance) { r.converged = true; break; }
  }
  r.beta = beta;
  r.log_partial_likelihood = st.loglik;
  r.covariance = detail::inverse(st.information);
  const double q = detail::z_two_sided(o.confidence);
  for (std::size_t j = 0; j < k; j++) {
    r.se.push_back(std::sqrt(r.covariance(j, j)));
    r.z.push_back(beta[j] / r.se[j]);
    r.p.push_back(2 * tyche::probability::Normal{0, 1}.survival(std::fabs(r.z[j])));
    r.ci_lower.push_back(beta[j] - q * r.se[j]);
    r.ci_upper.push_back(beta[j] + q * r.se[j]);
  }
  // Breslow's baseline cumulative hazard: Σ_{t_j ≤ t} d_j / Σ_{R(t_j)} exp(x β̂).
  std::map<double, std::pair<double, double>> by_time;  // time → (events, Σ exp η leaving at it)
  double s0 = 0;
  for (std::size_t i = 0; i < n; i++) {
    double eta = 0;
    for (std::size_t j = 0; j < k; j++) eta += X(i, j) * beta[j];
    by_time[time[i]].first += event[i] ? 1 : 0;
    by_time[time[i]].second += std::exp(eta);
    s0 += std::exp(eta);
  }
  double H = 0;
  for (const auto& [t, e] : by_time) {
    if (e.first > 0) { H += e.first / s0; r.baseline_time.push_back(t); r.baseline_hazard.push_back(H); }
    s0 -= e.second;
  }
  return r;
}

struct CumulativeIncidence {
  std::vector<double> time, incidence;
};

// cause: 0 censored, 1, 2, … the cause of the event; the incidence of `of_cause` at each time some event occurs.
inline CumulativeIncidence aalen_johansen(std::span<const double> time, std::span<const int> cause, int of_cause) {
  detail::require(time.size() == cause.size(), "one cause per time");
  std::map<double, std::pair<double, std::pair<double, double>>> by_time;  // time → (leaving, (all events, events of the cause))
  for (std::size_t i = 0; i < time.size(); i++) {
    auto& e = by_time[time[i]];
    e.first++;
    if (cause[i] != 0) e.second.first++;
    if (cause[i] == of_cause) e.second.second++;
  }
  CumulativeIncidence ci;
  double n = static_cast<double>(time.size()), S = 1, F = 0;
  for (const auto& [t, e] : by_time) {
    const double d = e.second.first, dc = e.second.second;
    if (d > 0) {
      F += S * dc / n;
      S *= 1 - d / n;
      ci.time.push_back(t);
      ci.incidence.push_back(F);
    }
    n -= e.first;
  }
  return ci;
}

}  // namespace peira::survival
