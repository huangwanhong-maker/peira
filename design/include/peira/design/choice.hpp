// Stated-choice designs: which alternatives to show together, so that the choices made among them say as much as possible about
// preferences (peira's decision 0007).
//
// THE DESIGN. Attributes a = 1 … A each take levels: a categorical attribute of L levels is coded by L − 1 columns of effects coding
// (the last level −1 in every column), a numeric one by its value. A profile is a level for each attribute; a choice set is J
// profiles; a design is S choice sets. Its coded matrix has one row of K columns per alternative, K the sum of the attributes' columns.
//
// ITS INFORMATION (Huber and Zwerina 1996, J. Marketing Research 33, 307–317; Kessels, Goos and Vandebroek 2006, J. Marketing
// Research 43, 409–419). Under the conditional logit with coefficients β, the information of a design is
//   M(β) = Σ_s Σ_j P_sj (x_sj − x̄_s)(x_sj − x̄_s)ᵀ,  x̄_s = Σ_j P_sj x_sj,
// and its D-error is det(M⁻¹)^{1/K}: smaller is better. With β = 0 it is the utility-neutral D₀-error (Street and Burgess 2007, The
// Construction of Optimal Stated Choice Experiments, Wiley); with priors β it is D_p; averaged over draws of a prior distribution,
// the Bayesian D_B-error (Sándor and Wedel 2001, J. Marketing Research 38, 430–444).
//
// A BOUND that the checks hold exactly. For binary attributes coded ±1, pairs and β = 0, a set's information is ddᵀ/4 with d = x₁ − x₂,
// so each diagonal entry is at most 1 per set and det M ≤ S^K: D₀ ≥ 1/S. The foldover of an orthogonal array (each set's second profile
// the first with every attribute changed) attains it.
//
// THE CONSTRUCTION. Coordinate exchange (Meyer and Nachtsheim 1995, Technometrics 37, 60–69; for choice designs, Kessels, Jones, Goos
// and Vandebroek 2009, J. Business Econ. Stat. 27, 279–291): from a random design, each attribute of each alternative of each set in
// turn is given the level that lowers the D-error most, until a pass changes nothing; the best of several random starts is kept. A set
// in which two alternatives are the same profile is never chosen (it says nothing). Draws come from the caller's calx stream.

#pragma once

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace peira::design {

struct ChoiceAttribute {
  std::string name;
  std::size_t levels = 2;                // categorical: L levels, coded by effects
  std::vector<double> values;            // numeric: the values of its levels (one column); empty for categorical
  std::size_t columns() const { return values.empty() ? levels - 1 : 1; }
  std::size_t count() const { return values.empty() ? levels : values.size(); }
};

struct ChoiceDesign {
  std::vector<ChoiceAttribute> attributes;
  std::size_t sets = 0, alternatives = 0;
  std::vector<std::size_t> level;        // [(s · J + j) · A + a]: the level of attribute a of alternative j of set s
  double d_error = std::numeric_limits<double>::infinity();
  std::size_t passes = 0;
  std::size_t at(std::size_t s, std::size_t j, std::size_t a) const { return level[(s * alternatives + j) * attributes.size() + a]; }
};

namespace detail_choice {

inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::design (choice): ") + what);
}

inline std::size_t width(const std::vector<ChoiceAttribute>& at) {
  std::size_t K = 0;
  for (const auto& a : at) K += a.columns();
  return K;
}

// The coded row of a profile.
inline void code(const std::vector<ChoiceAttribute>& at, const std::size_t* lv, double* row) {
  std::size_t c = 0;
  for (std::size_t a = 0; a < at.size(); a++) {
    const auto& A = at[a];
    if (!A.values.empty()) {
      row[c++] = A.values[lv[a]];
      continue;
    }
    for (std::size_t q = 0; q + 1 < A.levels; q++) row[c + q] = lv[a] == q ? 1.0 : (lv[a] == A.levels - 1 ? -1.0 : 0.0);
    c += A.levels - 1;
  }
}

}  // namespace detail_choice

// The information M(β) of a design.
inline calx::linalg::Matrix choice_information(const ChoiceDesign& d, const std::vector<double>& beta) {
  namespace dc = detail_choice;
  const std::size_t K = dc::width(d.attributes), J = d.alternatives, A = d.attributes.size();
  dc::require(beta.size() == K, "a prior for every coded column");
  calx::linalg::Matrix M(K);
  std::vector<double> X(J * K), P(J), mean(K);
  for (std::size_t s = 0; s < d.sets; s++) {
    double top = -std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < J; j++) {
      dc::code(d.attributes, &d.level[(s * J + j) * A], &X[j * K]);
      double v = 0;
      for (std::size_t c = 0; c < K; c++) v += X[j * K + c] * beta[c];
      P[j] = v, top = std::max(top, v);
    }
    double z = 0;
    for (double& p : P) p = std::exp(p - top), z += p;
    std::fill(mean.begin(), mean.end(), 0.0);
    for (std::size_t j = 0; j < J; j++) {
      P[j] /= z;
      for (std::size_t c = 0; c < K; c++) mean[c] += P[j] * X[j * K + c];
    }
    for (std::size_t j = 0; j < J; j++)
      for (std::size_t a = 0; a < K; a++)
        for (std::size_t b = 0; b < K; b++) M(a, b) += P[j] * (X[j * K + a] - mean[a]) * (X[j * K + b] - mean[b]);
  }
  return M;
}

// det(M⁻¹)^{1/K}; infinity where M is singular.
inline double d_error(const ChoiceDesign& d, const std::vector<double>& beta) {
  const auto M = choice_information(d, beta);
  const auto L = calx::linalg::log_determinant(M);
  if (L.sign <= 0 || !std::isfinite(L.log_abs)) return std::numeric_limits<double>::infinity();
  return std::exp(-L.log_abs / static_cast<double>(M.n));
}

// The mean D-error over draws of the prior (each draw a vector of K coefficients).
inline double bayesian_d_error(const ChoiceDesign& d, const std::vector<std::vector<double>>& draws) {
  double s = 0;
  for (const auto& b : draws) s += d_error(d, b);
  return s / static_cast<double>(draws.size());
}

// The design of `sets` sets of `alternatives` profiles by coordinate exchange from `starts` random starts, minimising the D-error at
// the prior β, or the Bayesian D-error over `draws` if they are given.
inline ChoiceDesign coordinate_exchange(const std::vector<ChoiceAttribute>& attributes, std::size_t sets, std::size_t alternatives, const std::vector<double>& beta,
                                        calx::random::Stream& stream, std::size_t starts = 10, const std::vector<std::vector<double>>& draws = {}) {
  namespace dc = detail_choice;
  dc::require(!attributes.empty() && sets >= 1 && alternatives >= 2, "attributes, a set and two alternatives at least");
  for (const auto& a : attributes) dc::require(a.count() >= 2, "two levels at least for every attribute");
  const std::size_t A = attributes.size(), J = alternatives;
  dc::require(dc::width(attributes) <= sets * (J - 1), "no more coded columns than the design can identify");
  const auto score = [&](const ChoiceDesign& d) {
    // A set with two identical profiles is never kept.
    for (std::size_t s = 0; s < d.sets; s++)
      for (std::size_t i = 0; i < J; i++)
        for (std::size_t j = i + 1; j < J; j++) {
          bool same = true;
          for (std::size_t a = 0; a < A && same; a++) same = d.at(s, i, a) == d.at(s, j, a);
          if (same) return std::numeric_limits<double>::infinity();
        }
    return draws.empty() ? d_error(d, beta) : bayesian_d_error(d, draws);
  };
  ChoiceDesign best;
  for (std::size_t r = 0; r < starts; r++) {
    ChoiceDesign d;
    d.attributes = attributes, d.sets = sets, d.alternatives = J;
    d.level.resize(sets * J * A);
    for (std::size_t s = 0; s < sets; s++)
      for (std::size_t j = 0; j < J; j++)
        for (std::size_t a = 0; a < A; a++) d.level[(s * J + j) * A + a] = static_cast<std::size_t>(stream.below(attributes[a].count()));
    double current = score(d);
    bool changed = true;
    while (changed && d.passes < 100) {
      changed = false;
      d.passes++;
      for (std::size_t cell = 0; cell < d.level.size(); cell++) {
        const std::size_t a = cell % A, keep = d.level[cell];
        std::size_t pick = keep;
        for (std::size_t v = 0; v < attributes[a].count(); v++) {
          if (v == keep) continue;
          d.level[cell] = v;
          const double e = score(d);
          if (e < current * (1 - 1e-12)) current = e, pick = v;
        }
        d.level[cell] = pick;
        changed = changed || pick != keep;
      }
    }
    d.d_error = current;
    if (d.d_error < best.d_error) best = d;
  }
  return best;
}

// The count of each level of each attribute over the design (level balance), [a][level].
inline std::vector<std::vector<std::size_t>> level_counts(const ChoiceDesign& d) {
  std::vector<std::vector<std::size_t>> out;
  for (const auto& a : d.attributes) out.emplace_back(a.count(), 0);
  for (std::size_t s = 0; s < d.sets; s++)
    for (std::size_t j = 0; j < d.alternatives; j++)
      for (std::size_t a = 0; a < d.attributes.size(); a++) out[a][d.at(s, j, a)]++;
  return out;
}

}  // namespace peira::design
