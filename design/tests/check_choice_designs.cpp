// Checks of stated-choice designs (choice.hpp).
//
// - The bound, exactly: three binary attributes in pairs, the foldover of the orthogonal array of four runs gives D₀ = 1/S; its
//   information is S I.
// - Coordinate exchange finds a design of 8 pairs at that bound, D₀ = 1/8, from random starts.
// - A set's information for pairs at β = 0 is ddᵀ/4 with d the difference of the two profiles.
// - With priors, the design made for them has a D_p-error no larger than the utility-neutral design's at those priors.
// - Level counts sum to the alternatives shown; a numeric attribute (a price of three levels) codes as one column.

#include <peira/design/choice.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace pd = peira::design;

int main() {
  calx::verification::Suite s;
  const std::vector<pd::ChoiceAttribute> binary{{"a", 2, {}}, {"b", 2, {}}, {"c", 2, {}}};
  {
    // The orthogonal array of four runs, 2^{3−1}, and its foldover: level 0 codes +1, level 1 codes −1.
    const std::size_t oa[4][3] = {{0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
    pd::ChoiceDesign d;
    d.attributes = binary, d.sets = 4, d.alternatives = 2;
    for (const auto& row : oa) {
      for (std::size_t a = 0; a < 3; a++) d.level.push_back(row[a]);
      for (std::size_t a = 0; a < 3; a++) d.level.push_back(1 - row[a]);
    }
    const std::vector<double> zero(3, 0.0);
    s.check("the foldover of the orthogonal array: D₀ = 1/S exactly", pd::d_error(d, zero) - 0.25, 1e-14);
    const auto M = pd::choice_information(d, zero);
    double off = 0;
    for (std::size_t a = 0; a < 3; a++)
      for (std::size_t b = 0; b < 3; b++) off = std::max(off, std::fabs(M(a, b) - (a == b ? 4.0 : 0.0)));
    s.check("... its information is S I", off, 1e-14);
  }
  {
    calx::random::Stream st(20261109, 0x63686f69u, 0, 0);  // "choi"
    const std::vector<double> zero(3, 0.0);
    const auto d = pd::coordinate_exchange(binary, 8, 2, zero, st, 20);
    s.check("coordinate exchange, 8 pairs of three binary attributes: the bound D₀ = 1/8 is reached", d.d_error - 0.125, 1e-12);
    // A set's information is ddᵀ/4.
    pd::ChoiceDesign one = d;
    one.sets = 1;
    one.level.resize(6);
    const auto M = pd::choice_information(one, zero);
    double worst = 0;
    for (std::size_t a = 0; a < 3; a++)
      for (std::size_t b = 0; b < 3; b++) {
        const double da = (one.level[a] == 0 ? 1.0 : -1.0) - (one.level[3 + a] == 0 ? 1.0 : -1.0);
        const double db = (one.level[b] == 0 ? 1.0 : -1.0) - (one.level[3 + b] == 0 ? 1.0 : -1.0);
        worst = std::max(worst, std::fabs(M(a, b) - da * db / 4));
      }
    s.check("a pair's information at β = 0 is ddᵀ/4", worst, 1e-14);
    const auto counts = pd::level_counts(d);
    bool sums = true;
    for (const auto& c : counts) sums = sums && c[0] + c[1] == 16;
    s.require("level counts sum to the 16 alternatives shown, for every attribute", sums);
  }
  {
    // A lane (binary), trees (three levels) and a price (numeric: 1, 2, 4), triples; priors (0.6; 0.3, 0.1; −0.5).
    const std::vector<pd::ChoiceAttribute> at{{"lane", 2, {}}, {"trees", 3, {}}, {"price", 3, {1, 2, 4}}};
    const std::vector<double> prior{0.6, 0.3, 0.1, -0.5}, zero(4, 0.0);
    calx::random::Stream st(20261109, 0x70726972u, 0, 0);  // "prir"
    const auto dp = pd::coordinate_exchange(at, 12, 3, prior, st, 10);
    const auto d0 = pd::coordinate_exchange(at, 12, 3, zero, st, 10);
    s.require("a numeric attribute codes as one column: four columns for a lane, trees and a price", pd::detail_choice::width(at) == 4);
    s.require("the design made for the priors is at least as good at them as the utility-neutral one (D_p " + std::to_string(pd::d_error(dp, prior)) + " against " +
                  std::to_string(pd::d_error(d0, prior)) + ")",
              pd::d_error(dp, prior) <= pd::d_error(d0, prior) * (1 + 1e-12));
    std::vector<std::vector<double>> draws;
    for (int r = 0; r < 20; r++) draws.push_back({0.6 + 0.2 * st.normal(), 0.3 + 0.1 * st.normal(), 0.1 + 0.1 * st.normal(), -0.5 + 0.2 * st.normal()});
    const double db = pd::bayesian_d_error(dp, draws);
    s.require("... and its Bayesian D-error over 20 draws of the prior is finite and positive: " + std::to_string(db), std::isfinite(db) && db > 0);
  }
  return s.done();
}
