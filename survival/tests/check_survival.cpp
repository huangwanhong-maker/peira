// Checks of peira/survival, against statsmodels 0.15.0 and lifelines 0.30.3 on Rossi's recidivism study and Walton's
// lifetimes (tests/fixtures/make_survival_fixtures.py writes the data and the values), and against exact identities.
//
// - Kaplan–Meier on Rossi: every estimate, Greenwood standard error and log(−log) interval at the 49 event times;
//   Walton's median, 56.
// - Nelson–Aalen on Rossi without tie smoothing.
// - The log-rank test of Walton's two groups (χ² 122.25) and of Rossi's financial aid (χ² 3.84, p 0.050).
// - Cox's model on Rossi's seven covariates with Breslow's and Efron's ties: coefficients, standard errors and the
//   partial log-likelihood, against statsmodels' PHReg; Efron's also against lifelines' CoxPHFitter.
// - Exact identities: with one covariate and no ties, Cox's score test at β = 0 is the log-rank statistic; Breslow's and
//   Efron's fits agree when no times are tied; the Aalen–Johansen incidences of every cause sum to 1 − Kaplan–Meier at
//   each time, and with one cause the incidence is 1 − Kaplan–Meier.

#include <peira/survival/survival.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace peira::survival;

namespace {

std::map<std::string, std::vector<double>> read_values(const std::string& path) {
  std::map<std::string, std::vector<double>> out;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string key;
    ss >> key;
    double v;
    while (ss >> v) out[key].push_back(v);
  }
  return out;
}

std::vector<std::vector<double>> read_csv(const std::string& path) {
  std::vector<std::vector<double>> rows;
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);  // the header
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::vector<double> row;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) row.push_back(std::stod(cell));
    rows.push_back(row);
  }
  return rows;
}

double worst_relative(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return INFINITY;
  double w = 0;
  for (std::size_t i = 0; i < a.size(); i++) w = std::max(w, std::fabs(a[i] - b[i]) / std::max(1e-300, std::fabs(b[i])));
  return w;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const std::string dir = PEIRA_FIXTURES;
  const auto V = read_values(dir + "/values.txt");
  const auto rossi = read_csv(dir + "/rossi.csv");
  const auto waltons = read_csv(dir + "/waltons.csv");
  if (V.empty() || rossi.size() != 432 || waltons.empty()) {
    s.require("fixtures: values.txt, rossi.csv (432 rows) and waltons.csv read", false, dir);
    return s.done();
  }
  std::vector<double> week, wt;
  std::vector<int> arrest, fin, we, wg;
  Rectangular X(432, 7);
  for (std::size_t i = 0; i < 432; i++) {
    week.push_back(rossi[i][0]);
    arrest.push_back(static_cast<int>(rossi[i][1]));
    fin.push_back(static_cast<int>(rossi[i][2]));
    for (int j = 0; j < 7; j++) X(i, j) = rossi[i][2 + j];
  }
  for (const auto& r : waltons) { wt.push_back(r[0]); we.push_back(static_cast<int>(r[1])); wg.push_back(static_cast<int>(r[2])); }

  {
    const KaplanMeier km = kaplan_meier(week, arrest);
    s.require("kaplan-meier (Rossi): the 49 event times", km.time == V.at("km_rossi.time"), std::to_string(km.time.size()));
    s.check("kaplan-meier (Rossi): S(t) against statsmodels (relative)", worst_relative(km.survival, V.at("km_rossi.survival")), 1e-14);
    s.check("kaplan-meier (Rossi): Greenwood's standard error against statsmodels (relative)", worst_relative(km.se, V.at("km_rossi.se")), 1e-12);
    s.check("kaplan-meier (Rossi): the log(−log) interval's lower bound against lifelines (relative)", worst_relative(km.lower, V.at("km_rossi.lower")), 1e-12);
    s.check("kaplan-meier (Rossi): its upper bound against lifelines (relative)", worst_relative(km.upper, V.at("km_rossi.upper")), 1e-12);
    const KaplanMeier kw = kaplan_meier(wt, we);
    s.check("kaplan-meier (Walton): the median against lifelines (56)", kw.median - V.at("km_waltons_median")[0], 0);
    const NelsonAalen na = nelson_aalen(week, arrest);
    s.check("nelson-aalen (Rossi): H(t) against lifelines (relative)", worst_relative(na.hazard, V.at("na_rossi.hazard")), 1e-13);
  }
  {
    const LogRank w = log_rank(wt, we, wg), r = log_rank(week, arrest, fin);
    s.check("log-rank (Walton, two groups): χ² against statsmodels (relative)", std::fabs(w.statistic / V.at("logrank_waltons.statistic")[0] - 1), 1e-12);
    s.check("log-rank (Walton): p against lifelines (relative)", std::fabs(w.p_value / V.at("logrank_waltons.lifelines_p_value")[0] - 1), 1e-8);
    s.check("log-rank (Rossi, financial aid): χ² against statsmodels (relative)", std::fabs(r.statistic / V.at("logrank_rossi_fin.statistic")[0] - 1), 1e-12);
    s.check("log-rank (Rossi): p against statsmodels (relative)", std::fabs(r.p_value / V.at("logrank_rossi_fin.p_value")[0] - 1), 1e-10);
  }
  for (const std::string ties : {"breslow", "efron"}) {
    CoxOptions o;
    o.ties = ties == "breslow" ? Ties::breslow : Ties::efron;
    const Cox c = cox(X, week, arrest, o);
    s.require("cox (Rossi, " + ties + "): converges", c.converged, std::to_string(c.iterations) + " iterations");
    s.check("cox (Rossi, " + ties + "): β against statsmodels (relative)", worst_relative(c.beta, V.at("cox_rossi_" + ties + ".beta")), 1e-8);
    s.check("cox (Rossi, " + ties + "): standard errors against statsmodels (relative)", worst_relative(c.se, V.at("cox_rossi_" + ties + ".se")), 1e-7);
    s.check("cox (Rossi, " + ties + "): partial log-likelihood against statsmodels", c.log_partial_likelihood - V.at("cox_rossi_" + ties + ".log_partial_likelihood")[0], 1e-8);
    if (ties == "efron") {
      s.check("cox (Rossi, efron): β against lifelines (relative)", worst_relative(c.beta, V.at("cox_rossi_efron_lifelines.beta")), 1e-7);
      s.check("cox (Rossi, efron): standard errors against lifelines (relative)", worst_relative(c.se, V.at("cox_rossi_efron_lifelines.se")), 1e-7);
    }
  }
  {
    // Untied times: Breslow = Efron, and the score test at 0 of one binary covariate is the log-rank test.
    calx::random::Stream st(41, 7, 0, 0);
    const std::size_t n = 120;
    std::vector<double> t(n);
    std::vector<int> e(n), g(n);
    Rectangular Z(n, 1);
    for (std::size_t i = 0; i < n; i++) {
      g[i] = static_cast<int>(st.below(2));
      t[i] = -std::log(st.uniform()) / (g[i] ? 1.6 : 1.0);
      const double c = -std::log(st.uniform()) / 0.4;
      e[i] = t[i] <= c;
      t[i] = std::min(t[i], c);
      Z(i, 0) = g[i];
    }
    CoxOptions ob, oe;
    ob.ties = Ties::breslow, oe.ties = Ties::efron;
    const Cox cb = cox(Z, t, e, ob), ce = cox(Z, t, e, oe);
    s.check("cox: with no tied times, Breslow's fit is Efron's", std::fabs(cb.beta[0] - ce.beta[0]), 1e-12);
    const auto st0 = detail::cox_state(Z, t, e, std::vector<double>{0.0}, Ties::breslow);
    const double score = st0.gradient[0] * st0.gradient[0] / st0.information(0, 0);
    std::vector<int> grp(g.begin(), g.end());
    for (auto& v : grp) v = 1 - v;  // the log-rank's first group, the treated
    s.check("cox: with no ties, the score test at β = 0 is the log-rank statistic", std::fabs(score - log_rank(t, e, grp).statistic), 1e-10);
  }
  {
    // Competing risks: three causes and censoring, with ties.
    calx::random::Stream st(43, 8, 0, 0);
    const std::size_t n = 300;
    std::vector<double> t(n);
    std::vector<int> cause(n), any(n);
    for (std::size_t i = 0; i < n; i++) {
      t[i] = std::round(-std::log(st.uniform()) * 10);  // ties on the whole numbers
      const double u = st.uniform();
      cause[i] = u < 0.15 ? 0 : u < 0.55 ? 1 : u < 0.85 ? 2 : 3;
      any[i] = cause[i] != 0;
    }
    const KaplanMeier km = kaplan_meier(t, any);
    const CumulativeIncidence f1 = aalen_johansen(t, cause, 1), f2 = aalen_johansen(t, cause, 2), f3 = aalen_johansen(t, cause, 3);
    double worst = 0;
    for (std::size_t j = 0; j < km.time.size(); j++) worst = std::max(worst, std::fabs(f1.incidence[j] + f2.incidence[j] + f3.incidence[j] - (1 - km.survival[j])));
    s.check("aalen-johansen: the causes' incidences sum to 1 − S at every time", worst, 1e-14);
    std::vector<int> one(n);
    for (std::size_t i = 0; i < n; i++) one[i] = any[i];
    const CumulativeIncidence fa = aalen_johansen(t, one, 1);
    double w1 = 0;
    for (std::size_t j = 0; j < km.time.size(); j++) w1 = std::max(w1, std::fabs(fa.incidence[j] - (1 - km.survival[j])));
    s.check("aalen-johansen: with one cause, the incidence is 1 − Kaplan–Meier", w1, 1e-14);
  }
  return s.done();
}
