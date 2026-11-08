// Designs: the assignment mechanisms by which units come to be treated.
//
// A design for n units is a law P(Z = z) on assignments z ∈ {0, 1}ⁿ (z_i = 1: unit i treated). Its assignment
// probabilities e_i = P(Z_i = 1) are what every estimator of peira divides by; a design is *known* when the law is the
// experimenter's, which is what makes a randomised experiment identify its effects without a model (docs/theory.md, §3).
//
// - Bernoulli: units treated independently, each with probability p.
// - Complete randomisation: exactly n₁ of n treated, every subset of size n₁ equally likely (1 / C(n, n₁)).
// - Blocked (stratified): complete randomisation within each block, independently across blocks.
// - Clustered: whole clusters assigned by complete randomisation of the clusters.
//
// Draws come from calx's counter-based streams, so a design drawn with a seed is the same everywhere. each_complete
// visits every assignment of complete randomisation, in lexicographic order, for exact randomisation distributions.
//
// Sources: Imbens and Rubin (2015), Causal Inference for Statistics, Social, and Biomedical Sciences, ch. 3–4 and 9;
// Fisher (1935), The Design of Experiments.

#pragma once

#include <calx/random/stream.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace peira::design {

using calx::random::Stream;
using Assignment = std::vector<int>;  // z_i ∈ {0, 1}

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::design: ") + what);
}
}  // namespace detail

// C(n, k) as a double (exact while it is below 2⁵³).
inline double binomial_coefficient(std::size_t n, std::size_t k) {
  if (k > n) return 0;
  if (k > n - k) k = n - k;
  double c = 1;
  for (std::size_t i = 1; i <= k; i++) c = c * static_cast<double>(n - k + i) / static_cast<double>(i);
  return std::round(c);
}

// Each unit treated independently with probability p.
inline Assignment bernoulli(std::size_t n, double p, Stream& s) {
  detail::require(p >= 0 && p <= 1, "a probability p in [0, 1]");
  Assignment z(n);
  for (auto& zi : z) zi = s.uniform() < p ? 1 : 0;
  return z;
}

// A uniformly random subset of size n1 of {0, …, n−1}, by a partial Fisher–Yates shuffle.
inline Assignment complete(std::size_t n, std::size_t n1, Stream& s) {
  detail::require(n1 <= n, "at most n of n units treated");
  std::vector<std::size_t> index(n);
  for (std::size_t i = 0; i < n; i++) index[i] = i;
  for (std::size_t i = 0; i < n1; i++) {
    const std::size_t j = i + static_cast<std::size_t>(s.below(n - i));
    std::swap(index[i], index[j]);
  }
  Assignment z(n, 0);
  for (std::size_t i = 0; i < n1; i++) z[index[i]] = 1;
  return z;
}

// Complete randomisation within each block: treated[b] of the units labelled b are treated.
inline Assignment blocked(std::span<const int> block, const std::map<int, std::size_t>& treated, Stream& s) {
  std::map<int, std::vector<std::size_t>> members;
  for (std::size_t i = 0; i < block.size(); i++) members[block[i]].push_back(i);
  Assignment z(block.size(), 0);
  for (const auto& [b, units] : members) {
    const auto it = treated.find(b);
    detail::require(it != treated.end(), "a number treated for every block");
    const Assignment zb = complete(units.size(), it->second, s);
    for (std::size_t k = 0; k < units.size(); k++) z[units[k]] = zb[k];
  }
  return z;
}

// Whole clusters treated: clusters_treated of the distinct cluster labels, by complete randomisation of the clusters.
inline Assignment clustered(std::span<const int> cluster, std::size_t clusters_treated, Stream& s) {
  std::map<int, std::size_t> label;
  for (int c : cluster) label.emplace(c, 0);
  std::size_t k = 0;
  for (auto& [c, index] : label) index = k++;
  const Assignment zc = complete(label.size(), clusters_treated, s);
  Assignment z(cluster.size());
  for (std::size_t i = 0; i < cluster.size(); i++) z[i] = zc[label[cluster[i]]];
  return z;
}

// Every assignment of complete randomisation (n1 of n treated), in lexicographic order of the treated sets; f(z) for each.
// There are C(n, n1) of them: the caller bounds n.
template <class F>
void each_complete(std::size_t n, std::size_t n1, F&& f) {
  detail::require(n1 <= n, "at most n of n units treated");
  std::vector<std::size_t> c(n1);
  for (std::size_t i = 0; i < n1; i++) c[i] = i;
  Assignment z(n);
  while (true) {
    std::fill(z.begin(), z.end(), 0);
    for (std::size_t i : c) z[i] = 1;
    f(static_cast<const Assignment&>(z));
    // The next combination.
    std::size_t i = n1;
    while (i > 0 && c[i - 1] == n - n1 + i - 1) i--;
    if (i == 0) return;
    c[i - 1]++;
    for (std::size_t j = i; j < n1; j++) c[j] = c[j - 1] + 1;
  }
}

// Every assignment of a Bernoulli design with unit-specific probabilities e_i, with its probability: f(z, P(z)).
// There are 2ⁿ: the caller bounds n.
template <class F>
void each_bernoulli(std::span<const double> e, F&& f) {
  const std::size_t n = e.size();
  detail::require(n < 31, "fewer than 31 units to enumerate");
  Assignment z(n);
  for (std::uint32_t mask = 0; mask < (1u << n); mask++) {
    double p = 1;
    for (std::size_t i = 0; i < n; i++) {
      z[i] = (mask >> i) & 1u;
      p *= z[i] ? e[i] : 1 - e[i];
    }
    f(static_cast<const Assignment&>(z), p);
  }
}

}  // namespace peira::design
