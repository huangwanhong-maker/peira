// Checks of peira/detection: distance sampling (distance.hpp).
//
// - Closed forms: without truncation the half-normal's MLE is σ̂² = Σx²/n (lines) and Σr²/(2n) (points), its effective strip
//   half-width σ√(π/2) and effective area 2πσ²; the delta method's cv of μ̂ is 1/√(2n) (lines) and of ν̂ 1/√n (points), from
//   the information of log σ.
// - The closed forms with truncation against calx's quadrature of g; the hazard rate's quadrature on [0, w] plus its tail
//   σ(w/σ)^{1−b}/(b − 1) against ∫₀^∞ g = σ Γ(1 − 1/b).
// - Horvitz–Thompson: D̂ times the covered area is peira::weighting's Horvitz–Thompson total of the detections with the
//   estimated P_a.
// - The encounter rate's variance (Fewster et al. 2009, R2) by hand on three lines.
// - Cosine adjustments keep g(0) = 1; a g negative within [0, w] is refused; on half-normal distances the adjustment of order 2
//   is within three standard errors of 0.
// - Recovery from simulated surveys: lines with a half-normal (200 surveys) and points with a hazard rate (150 surveys):
//   D̂'s bias below a quarter (lines) or 0.3 (points) of its spread, the log-normal intervals' coverage within 3.5 Monte Carlo
//   standard errors of 95%, the standard errors' mean within 15% (lines) or 20% (points) of the spread; AIC chooses the hazard
//   rate on a large hazard-rate sample.

#include <peira/detection/detection.hpp>
#include <peira/weighting/weighting.hpp>

#include <calx/quadrature/gauss_kronrod.hpp>
#include <calx/random/stream.hpp>
#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

namespace pd = peira::detection;

namespace {

struct Recovery {
  double sum = 0, sum2 = 0, se = 0, cover = 0, n = 0;
  void add(double d, double se_, double lo, double hi, double truth) {
    sum += d, sum2 += d * d, se += se_, cover += lo <= truth && truth <= hi, n++;
  }
  double mean() const { return sum / n; }
  double sd() const { return std::sqrt(std::max(0.0, sum2 / n - mean() * mean())); }
};

}  // namespace

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;

  // ---- The half-normal without truncation.
  {
    calx::random::Stream st(2001, 1, 0, 0);
    const double sigma = 20;
    auto g = [&](double x) { return std::exp(-x * x / (2 * sigma * sigma)); };
    auto lines = pd::simulate_distances(st, pd::Transect::line, 0.002, std::vector<double>(10, 300.0), 12 * sigma, g);
    lines.truncation = std::numeric_limits<double>::infinity();
    const auto fl = pd::fit_detection(lines, {pd::Key::half_normal, {}});
    double sx2 = 0;
    for (double x : lines.distance) sx2 += x * x;
    const double n = static_cast<double>(lines.distance.size());
    s.check("half-normal lines, untruncated: σ̂² = Σx²/n (relative)", fl.sigma * fl.sigma / (sx2 / n) - 1, 1e-8);
    s.check("half-normal lines, untruncated: μ̂ = σ̂√(π/2) (relative)", fl.effective / (fl.sigma * std::sqrt(pi / 2)) - 1, 1e-14);
    s.check("half-normal lines, untruncated: cv(μ̂) = 1/√(2n) by the delta method (relative)", fl.effective_cv * std::sqrt(2 * n) - 1, 1e-5);
    auto points = pd::simulate_distances(st, pd::Transect::point, 0.002, std::vector<double>(40, 1.0), 12 * sigma, g);
    points.truncation = std::numeric_limits<double>::infinity();
    const auto fp = pd::fit_detection(points, {pd::Key::half_normal, {}});
    double sr2 = 0;
    for (double r : points.distance) sr2 += r * r;
    const double m = static_cast<double>(points.distance.size());
    s.check("half-normal points, untruncated: σ̂² = Σr²/(2n) (relative)", fp.sigma * fp.sigma / (sr2 / (2 * m)) - 1, 1e-8);
    s.check("half-normal points, untruncated: ν̂ = 2πσ̂² (relative)", fp.effective / (2 * pi * fp.sigma * fp.sigma) - 1, 1e-14);
    s.check("half-normal points, untruncated: cv(ν̂) = 1/√n by the delta method (relative)", fp.effective_cv * std::sqrt(m) - 1, 1e-5);
  }
  // ---- Closed forms against quadrature; the hazard rate's tail.
  {
    const double sigma = 13, w = 30;
    const pd::DetectionModel hn{pd::Key::half_normal, {}};
    const std::vector<double> th{std::log(sigma)};
    auto g = [&](double x) { return std::exp(-x * x / (2 * sigma * sigma)); };
    const auto ql = calx::quadrature::integrate(g, 0.0, w, 1e-300, 1e-14);
    const auto qp = calx::quadrature::integrate([&](double r) { return r * g(r); }, 0.0, w, 1e-300, 1e-14);
    s.check("half-normal, truncated: σ√(π/2) erf(w/(σ√2)) against quadrature of g (relative)", pd::detail::effective(hn, pd::Transect::line, th, w) / ql.value - 1, 1e-13);
    s.check("half-normal, truncated: 2πσ²(1 − e^{−w²/(2σ²)}) against 2π ∫ r g (relative)", pd::detail::effective(hn, pd::Transect::point, th, w) / (2 * pi * qp.value) - 1, 1e-13);
    s.check("half-normal, untruncated: μ = σ√(π/2) against quadrature on [0, 40σ] (relative)",
            pd::detail::effective(hn, pd::Transect::line, th, std::numeric_limits<double>::infinity()) / calx::quadrature::integrate(g, 0.0, 40 * sigma, 1e-300, 1e-14).value - 1, 1e-13);
    const double b = 3, W = 50 * sigma;
    const pd::DetectionModel hr{pd::Key::hazard_rate, {}};
    const double mu = pd::detail::effective(hr, pd::Transect::line, std::vector<double>{std::log(sigma), std::log(b)}, W);
    const double tail = sigma * std::pow(W / sigma, 1 - b) / (b - 1);
    s.check("hazard rate: ∫₀^w g + σ(w/σ)^{1−b}/(b − 1) = σ Γ(1 − 1/b) (relative)", (mu + tail) / (sigma * std::exp(calx::special_functions::log_gamma(1 - 1 / b))) - 1, 1e-9);
  }
  // ---- Horvitz–Thompson, and the encounter rate by hand.
  {
    calx::random::Stream st(2001, 2, 0, 0);
    const double sigma = 15, w = 40;
    auto g = [&](double x) { return std::exp(-x * x / (2 * sigma * sigma)); };
    const auto survey = pd::simulate_distances(st, pd::Transect::line, 0.001, std::vector<double>(12, 400.0), w, g);
    const auto fit = pd::fit_detection(survey, {pd::Key::half_normal, {}});
    const auto d = pd::density(survey, fit);
    std::vector<double> ones(d.n, 1.0), pa(d.n, fit.probability);
    const double ht = peira::weighting::horvitz_thompson_total(ones, pa);
    s.check("horvitz-thompson: D̂ × 2wL = Σ 1/P̂_a over the detections (relative)", d.density * d.covered_area / ht - 1, 1e-13);
    s.check("horvitz-thompson: P̂_a = μ̂/w", fit.probability - fit.effective / w, 1e-15);

    pd::DistanceSurvey hand;
    hand.transect = pd::Transect::line, hand.truncation = 10, hand.effort = {2, 3, 5};
    const std::size_t per[3] = {4, 3, 9};
    for (std::size_t k = 0; k < 3; k++)
      for (std::size_t a = 0; a < per[k]; a++) hand.distance.push_back(0.5 + a), hand.on.push_back(k);
    const auto uf = pd::fit_detection(hand, {pd::Key::uniform, {}});
    const auto ud = pd::density(hand, uf);
    s.check("encounter rate: n/E = 16/10", ud.encounter_rate - 1.6, 1e-15);
    s.check("encounter rate: R2's variance K/(E²(K − 1)) Σ e_k²(n_k/e_k − n/E)² = 0.0732 by hand", ud.encounter_rate_se * ud.encounter_rate_se - 0.0732, 1e-15);
    s.check("uniform key: D̂ = (n/E)/(2w), with no variance from detection", std::fabs(ud.density - 1.6 / 20) + ud.cv_detection, 1e-15);
    s.check("uniform key: Satterthwaite's degrees of freedom are K − 1 = 2", ud.df - 2, 1e-12);
  }
  // ---- Cosine adjustments.
  {
    calx::random::Stream st(2001, 3, 0, 0);
    const double sigma = 15, w = 40;
    auto g = [&](double x) { return std::exp(-x * x / (2 * sigma * sigma)); };
    const auto survey = pd::simulate_distances(st, pd::Transect::line, 0.004, std::vector<double>(12, 400.0), w, g);
    const pd::DetectionModel adj{pd::Key::half_normal, {2}};
    const auto fit = pd::fit_detection(survey, adj);
    s.check("cosine adjustments: g(0) = 1", pd::detection_probability(fit, 0) - 1, 1e-15);
    s.require("cosine adjustments: converges, with a positive definite information", fit.converged && fit.information_positive);
    s.check("cosine adjustments: on half-normal distances a₂/se(a₂) (limit 3)", fit.theta[1] / fit.se[1], 3);
    const pd::DetectionModel uni{pd::Key::uniform, {1, 2}}, one{pd::Key::uniform, {1}};
    // g = (1 + 1.5 cos(πx/w))/2.5: positive at 0 and at every distance below 0.73w, negative at w.
    std::vector<double> near;
    for (double x : survey.distance)
      if (x < 0.7 * w) near.push_back(x);
    s.require("cosine adjustments: a g negative within [0, w] is refused, though positive at every distance seen",
              std::isinf(pd::detail::distance_log_likelihood(one, pd::Transect::line, near, w, std::vector<double>{1.5})) &&
                  std::isfinite(pd::detail::distance_log_likelihood(one, pd::Transect::line, near, w, std::vector<double>{0.5})));
    const auto fu = pd::fit_detection(survey, uni);
    const auto du = pd::density(survey, fu);
    const auto dh = pd::density(survey, pd::fit_detection(survey, {pd::Key::half_normal, {}}));
    s.check("cosine adjustments: the uniform key with cos(πx/w), cos(2πx/w) gives D̂ within 3 s.e. of the truth", (du.density - 0.004) / du.se, 3);
    s.check("cosine adjustments: … and within 2 s.e. of the half-normal's D̂", (du.density - dh.density) / du.se, 2);
  }
  // ---- Recovery: lines with a half-normal.
  {
    const double sigma = 15, w = 45, D = 0.0005;
    auto g = [&](double x) { return std::exp(-x * x / (2 * sigma * sigma)); };
    Recovery r;
    double worst_converged = 1;
    for (std::uint32_t rep = 0; rep < 200; rep++) {
      calx::random::Stream st(2001, 4, rep, 0);
      const auto survey = pd::simulate_distances(st, pd::Transect::line, D, std::vector<double>(20, 500.0), w, g);
      const auto fit = pd::fit_detection(survey, {pd::Key::half_normal, {}});
      if (!fit.converged) worst_converged = 0;
      const auto d = pd::density(survey, fit);
      r.add(d.density, d.se, d.lower, d.upper, D);
    }
    s.require("line recovery: every fit converges", worst_converged == 1);
    s.check("line recovery: D̂'s bias / its spread (Monte Carlo error 0.071)", (r.mean() - D) / r.sd(), 0.25);
    s.check("line recovery: the log-normal interval's coverage − 0.95 (limit 3.5 Monte Carlo s.e.)", r.cover / r.n - 0.95, 3.5 * std::sqrt(0.95 * 0.05 / r.n));
    s.check("line recovery: the mean standard error / D̂'s spread − 1", r.se / r.n / r.sd() - 1, 0.15);
  }
  // ---- Recovery: points with a hazard rate, visited twice.
  {
    const double sigma = 10, b = 3, w = 30, D = 0.005;
    auto g = [&](double x) { return -std::expm1(-std::pow(x / sigma, -b)); };
    Recovery r;
    for (std::uint32_t rep = 0; rep < 150; rep++) {
      calx::random::Stream st(2001, 5, rep, 0);
      const auto survey = pd::simulate_distances(st, pd::Transect::point, D, std::vector<double>(30, 2.0), w, g);
      const auto d = pd::density(survey, pd::fit_detection(survey, {pd::Key::hazard_rate, {}}));
      r.add(d.density, d.se, d.lower, d.upper, D);
    }
    s.check("point recovery: D̂'s bias / its spread (Monte Carlo error 0.082)", (r.mean() - D) / r.sd(), 0.3);
    s.check("point recovery: the log-normal interval's coverage − 0.95 (limit 3.5 Monte Carlo s.e.)", r.cover / r.n - 0.95, 3.5 * std::sqrt(0.95 * 0.05 / r.n));
    s.check("point recovery: the mean standard error / D̂'s spread − 1", r.se / r.n / r.sd() - 1, 0.2);
  }
  // ---- AIC between keys.
  {
    calx::random::Stream st(2001, 6, 0, 0);
    const double sigma = 10, b = 2.5, w = 40;
    auto g = [&](double x) { return -std::expm1(-std::pow(x / sigma, -b)); };
    const auto survey = pd::simulate_distances(st, pd::Transect::line, 0.01, std::vector<double>(10, 1000.0), w, g);
    const std::vector<pd::DetectionModel> models{{pd::Key::half_normal, {}}, {pd::Key::hazard_rate, {}}, {pd::Key::uniform, {1, 2}}};
    const auto a = pd::distance_sampling(survey, models);
    s.require("AIC: chooses the hazard rate on " + std::to_string(a.density.n) + " hazard-rate distances", a.chosen == 1,
              "ΔAIC " + calx::verification::exponential(a.delta_aic[0]) + ", " + calx::verification::exponential(a.delta_aic[2]));
    s.check("AIC: the chosen model's shape b̂ within three standard errors of 2.5 (on log b)", (a.fits[1].theta[1] - std::log(b)) / a.fits[1].se[1], 3);
    // One survey: a z-score, whose calibration over many surveys is the point recovery's above (on 60 surveys of this design
    // its mean was −0.04 and its spread 1.01).
    s.check("AIC: the chosen model's (D̂ − D)/se on this one survey (limit 4)", (a.density.density - 0.01) / a.density.se, 4);
  }
  return s.done();
}
