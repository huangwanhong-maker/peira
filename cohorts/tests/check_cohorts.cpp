// Checks of cohorts (cohorts.hpp).
//
// - One life on the Lexis diagram, split by hand: born 1950.5, followed from 2000 to 2003.25, it gives half-years to six cells and a
//   quarter to the seventh, and its death to the cell of age 52 in 2003.
// - Random lives: the person-time inside and outside the grid sums to their follow-up exactly, and every event is counted once; the
//   split agrees with a brute-force integration on a fine grid of instants.
// - Garwood's interval: at its bounds the Poisson's tails at D are α/2 exactly, for D = 0, 1, 5, 30.
// - Standardisation: equal rates in every stratum standardise to that rate, whatever the weights; observed events equal to those
//   expected give an SMR of one, inside its interval.

#include <peira/cohorts/cohorts.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace pc = peira::cohorts;

int main() {
  calx::verification::Suite s;
  pc::Grid g;
  for (int a = 49; a <= 53; a++) g.age.push_back(a);
  for (int p = 2000; p <= 2004; p++) g.period.push_back(p);
  {
    const auto L = pc::lexis({{1950.5, 2000.0, 2003.25, true}}, g);
    // [age band][period band]: ages 49..52, periods 2000..2003.
    const double expect[4][4] = {{0.5, 0, 0, 0}, {0.5, 0.5, 0, 0}, {0, 0.5, 0.5, 0}, {0, 0, 0.5, 0.25}};
    double worst = 0;
    for (int a = 0; a < 4; a++)
      for (int p = 0; p < 4; p++) worst = std::max(worst, std::fabs(L.person_time[a][p] - expect[a][p]));
    s.check("Lexis: one life born 1950.5, followed 2000–2003.25, split by hand into six half-years and a quarter", worst, 1e-12);
    s.require("... its death in the cell of age 52 in 2003", L.events[3][3] == 1 && L.events_outside == 0);
  }
  {
    calx::random::Stream st(20261108, 0x6c657869u, 0, 0);  // "lexi"
    std::vector<pc::Life> lives;
    double follow = 0, events = 0;
    for (int k = 0; k < 400; k++) {
      const double b = st.uniform(1940, 1960), e = st.uniform(1995, 2004), x = e + st.uniform(0, 8);
      const bool ev = st.uniform() < 0.3;
      lives.push_back({b, e, x, ev});
      follow += x - e, events += ev;
    }
    const auto L = pc::lexis(lives, g);
    double inside = 0, counted = 0;
    for (const auto& row : L.person_time)
      for (double v : row) inside += v;
    for (const auto& row : L.events)
      for (double v : row) counted += v;
    s.check("Lexis: 400 random lives, person-time inside and outside the grid sums to their follow-up (relative)", (inside + L.outside - follow) / follow, 1e-12);
    s.require("... and every event is counted once", counted + L.events_outside == events, std::to_string(counted) + " + " + std::to_string(L.events_outside));
    // Brute force: the instants of a fine grid, each assigned to its cell.
    std::vector<std::vector<double>> brute(4, std::vector<double>(4, 0.0));
    const double dt = 1e-3;
    for (const auto& life : lives)
      for (double t = life.entry + dt / 2; t < life.exit; t += dt) {
        const double age = t - life.birth;
        if (age < 49 || age >= 53 || t < 2000 || t >= 2004) continue;
        brute[static_cast<int>(age) - 49][static_cast<int>(t) - 2000] += dt;
      }
    double worst = 0;
    for (int a = 0; a < 4; a++)
      for (int p = 0; p < 4; p++) worst = std::max(worst, std::fabs(L.person_time[a][p] - brute[a][p]));
    s.check("... the split agrees with a brute-force integration over instants a thousandth of a year apart (years)", worst, 0.02);
  }
  {
    double worst = 0;
    for (const double d : {0.0, 1.0, 5.0, 30.0}) {
      const auto r = pc::rate(d, 10, 0.95);
      const double upper_tail = tyche::probability::Poisson{r.upper * 10}.cdf(static_cast<std::int64_t>(d));
      worst = std::max(worst, std::fabs(upper_tail - 0.025));
      if (d > 0) worst = std::max(worst, std::fabs(tyche::probability::Poisson{r.lower * 10}.survival(static_cast<std::int64_t>(d) - 1) - 0.025));
      else worst = std::max(worst, r.lower);
    }
    s.check("Garwood: at the bounds the Poisson's tails at D are 2.5% (D = 0, 1, 5, 30)", worst, 1e-9);
  }
  {
    const auto d = pc::direct({2, 6, 12}, {100, 300, 600}, {5, 1, 0.25});
    s.check("direct standardisation: equal rates (0.02) standardise to that rate, whatever the weights", d.value - 0.02, 1e-15);
    const auto i = pc::indirect(20, {100, 300, 600}, {0.01, 0.02, 13.0 / 600});
    s.check("indirect: observed events equal to those expected give an SMR of one", i.value - 1, 1e-12);
    s.require("... inside its interval", i.lower < 1 && i.upper > 1);
  }
  return s.done();
}
