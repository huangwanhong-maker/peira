// Checks of instrumental variables (instrumental.hpp), on simulated data with an endogenous
// regressor.
//
// Just identified, 2SLS is the simple IV estimator (ZᵀX)⁻¹Zᵀy and its classical covariance
// σ̂² (ZᵀX)⁻¹ ZᵀZ (XᵀZ)⁻¹. Overidentified with two valid instruments: OLS is biased, 2SLS recovers
// the coefficients within four standard errors; the first-stage F is large and the partial R²
// in (0, 1); Sargan's and Hansen's statistics do not reject. One instrument invalid (entering the
// outcome): both reject. A nearly irrelevant instrument: the first-stage F is small (below 10,
// Staiger and Stock 1997). Clustered and HC1 covariances are positive and differ.

#include <peira/quasi/instrumental.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace ne = peira::quasi;
namespace tr = tyche::regression;
using calx::linalg::Matrix;
using calx::linalg::Rectangular;
using calx::verification::exponential;
using calx::verification::relative;

namespace {

calx::random::Stream z(20261025, 0x69767276u, 0, 0);  // "ivrv"

struct Sample {
  std::vector<double> y;
  Rectangular W, E, Z;  // [1, w], [x], [z1, z2]
  std::vector<std::int64_t> group;
};

// y = 1 − w + 2x + u, x = π₁ z1 + π₂ z2 + 0.3 w + v, u = 0.6 v + e; `invalid` puts z2 into y.
Sample simulate(std::size_t n, double pi1, double pi2, double invalid) {
  Sample s;
  s.W = Rectangular(n, 2);
  s.E = Rectangular(n, 1);
  s.Z = Rectangular(n, 2);
  for (std::size_t i = 0; i < n; i++) {
    const double w = z.normal(), z1 = z.normal(), z2 = z.normal(), v = z.normal(), e = z.normal();
    const double x = pi1 * z1 + pi2 * z2 + 0.3 * w + v;
    const double u = 0.6 * v + e;
    s.y.push_back(1 - w + 2 * x + invalid * z2 + u);
    s.W(i, 0) = 1;
    s.W(i, 1) = w;
    s.E(i, 0) = x;
    s.Z(i, 0) = z1;
    s.Z(i, 1) = z2;
    s.group.push_back(static_cast<std::int64_t>(i % 50));
  }
  return s;
}

Rectangular first_column(const Rectangular& A) {
  Rectangular B(A.rows, 1);
  for (std::size_t i = 0; i < A.rows; i++) B(i, 0) = A(i, 0);
  return B;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // Just identified: the simple IV estimator.
  {
    const auto d = simulate(800, 0.8, 0.5, 0);
    const auto Z1 = first_column(d.Z);
    ne::IvOptions o;
    o.covariance = ne::IvCovariance::classical;
    const auto iv = ne::two_stage_least_squares(d.y, d.W, d.E, Z1, o);
    // (ZᵀX)⁻¹ Zᵀy with Z = [1, w, z1], X = [1, w, x].
    const std::size_t n = d.y.size();
    Matrix ZX(3), ZZ(3), Zy(3);
    for (std::size_t i = 0; i < n; i++) {
      const double zr[3] = {1, d.W(i, 1), d.Z(i, 0)}, xr[3] = {1, d.W(i, 1), d.E(i, 0)};
      for (int a = 0; a < 3; a++) {
        Zy(a, 0) += zr[a] * d.y[i];
        for (int b = 0; b < 3; b++) ZX(a, b) += zr[a] * xr[b], ZZ(a, b) += zr[a] * zr[b];
      }
    }
    const Matrix b = calx::linalg::solve(ZX, Zy);
    double db = 0;
    for (int a = 0; a < 3; a++) db = std::max(db, std::fabs(iv.beta[static_cast<std::size_t>(a)] - b(a, 0)));
    s.check("just identified: 2SLS is (ZᵀX)⁻¹Zᵀy", db, 1e-10);
    // σ̂² (ZᵀX)⁻¹ ZᵀZ (XᵀZ)⁻¹.
    double uu = 0;
    for (double u : iv.residuals) uu += u * u;
    Matrix XZ(3);
    for (int a = 0; a < 3; a++)
      for (int c = 0; c < 3; c++) XZ(a, c) = ZX(c, a);
    const Matrix inv = calx::linalg::solve(ZX, Matrix::identity(3)), invT = calx::linalg::solve(XZ, Matrix::identity(3));
    const Matrix V = (uu / (static_cast<double>(n) - 3)) * (inv * ZZ * invT);
    double dv = 0;
    for (int a = 0; a < 3; a++)
      for (int c = 0; c < 3; c++) dv = std::max(dv, relative(iv.covariance(a, c), V(a, c)));
    s.check("... and its classical covariance σ̂² (ZᵀX)⁻¹ZᵀZ(XᵀZ)⁻¹ (relative)", dv, 1e-8);
    s.require("... with no overidentifying restriction to test (NaN)", std::isnan(iv.sargan) && std::isnan(iv.hansen_j) && iv.overidentification_df == 0);
  }
  // Overidentified, valid instruments.
  {
    const auto d = simulate(3000, 0.8, 0.5, 0);
    const auto iv = ne::two_stage_least_squares(d.y, d.W, d.E, d.Z);
    const auto ols = tr::ols(Rectangular::from_columns({d.W.column(0), d.W.column(1), d.E.column(0)}), d.y, {});
    s.require("two valid instruments, n = 3000: OLS is biased (its estimate of 2 off by more than ten standard errors)", std::fabs(ols.beta[2] - 2) > 10 * ols.se[2],
              exponential(ols.beta[2]));
    const double t0 = std::fabs(iv.beta[0] - 1) / iv.se[0], t1 = std::fabs(iv.beta[1] + 1) / iv.se[1], t2 = std::fabs(iv.beta[2] - 2) / iv.se[2];
    s.require("... 2SLS recovers (1, −1, 2) within four standard errors", t0 < 4 && t1 < 4 && t2 < 4,
              exponential(iv.beta[2]) + " ± " + exponential(iv.se[2]));
    s.require("... the first-stage F is large and the partial R² in (0, 1)",
              iv.first_stage[0].f_statistic > 100 && iv.first_stage[0].partial_r_squared > 0 && iv.first_stage[0].partial_r_squared < 1,
              "F = " + exponential(iv.first_stage[0].f_statistic) + ", partial R² = " + exponential(iv.first_stage[0].partial_r_squared));
    s.require("... Sargan and Hansen do not reject (p > 0.01)", iv.sargan_p_value > 0.01 && iv.hansen_p_value > 0.01 && iv.overidentification_df == 1,
              "Sargan " + exponential(iv.sargan) + ", J " + exponential(iv.hansen_j));
    s.check("... and under homoskedasticity they nearly agree", std::fabs(iv.sargan - iv.hansen_j), 0.5);
    ne::IvOptions oc;
    oc.covariance = ne::IvCovariance::cluster;
    oc.clusters = d.group;
    const auto ivc = ne::two_stage_least_squares(d.y, d.W, d.E, d.Z, oc);
    s.require("... clustered in 50 groups: the same estimate, t on 49 degrees of freedom, positive variances",
              std::fabs(ivc.beta[2] - iv.beta[2]) < 1e-12 && ivc.df_inference == 49 && ivc.se[2] > 0 && ivc.se[2] != iv.se[2]);
  }
  // One invalid instrument.
  {
    const auto d = simulate(3000, 0.8, 0.5, 0.5);
    const auto iv = ne::two_stage_least_squares(d.y, d.W, d.E, d.Z);
    s.require("an invalid instrument (z₂ in the outcome): Sargan and Hansen reject (p < 10⁻⁶)", iv.sargan_p_value < 1e-6 && iv.hansen_p_value < 1e-6,
              "Sargan " + exponential(iv.sargan) + ", J " + exponential(iv.hansen_j));
  }
  // A weak instrument.
  {
    const auto d = simulate(2000, 0.02, 0.0, 0);
    const auto iv = ne::two_stage_least_squares(d.y, d.W, d.E, first_column(d.Z));
    s.require("a nearly irrelevant instrument: the first-stage F is below 10", iv.first_stage[0].f_statistic < 10, "F = " + exponential(iv.first_stage[0].f_statistic));
  }
  return s.done();
}
