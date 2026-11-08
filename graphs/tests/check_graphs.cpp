// Checks of causal graphs (causal.hpp).
//
// - The three elementary graphs: a chain and a fork are blocked by their middle; a collider is blocked by nothing and opened by its
//   middle or a descendant of it.
// - M-bias: X ← U₁ → M ← U₂ → Y with X → Y: the empty set satisfies the back-door criterion and {M} does not, {M, U₁} does; the
//   minimal sets among U₁, U₂ and M, by enumeration, are the empty set alone (every other valid set contains it).
// - The front door: U → X, U → Y, X → M → Y with U unobserved: {M} satisfies the front-door criterion, and no observed set the
//   back-door one.
// - Random graphs against linear Gaussian models with random coefficients (faithful almost surely): for every pair and every
//   conditioning set of at most two, d-separation holds exactly when the partial correlation of the model's covariance vanishes.
// - Acyclicity is recognised.

#include <peira/graphs/causal.hpp>

#include <calx/linalg/rectangular.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <set>
#include <string>
#include <vector>

namespace pg = peira::graphs;

namespace {
// The partial correlation of i and j given zs from a covariance matrix: from the inverse of the submatrix on {i, j} ∪ zs.
double partial_correlation(const std::vector<std::vector<double>>& S, std::size_t i, std::size_t j, const std::set<std::size_t>& zs) {
  std::vector<std::size_t> idx{i, j};
  idx.insert(idx.end(), zs.begin(), zs.end());
  const std::size_t m = idx.size();
  calx::linalg::Matrix A(m), I(m);
  for (std::size_t a = 0; a < m; a++) {
    I(a, a) = 1;
    for (std::size_t b = 0; b < m; b++) A(a, b) = S[idx[a]][idx[b]];
  }
  const auto P = calx::linalg::solve(A, I);
  return -P(0, 1) / std::sqrt(P(0, 0) * P(1, 1));
}
}  // namespace

int main() {
  calx::verification::Suite s;
  {
    pg::Dag chain(3), fork(3), collider(4);
    chain.arrow(0, 1), chain.arrow(1, 2);
    fork.arrow(1, 0), fork.arrow(1, 2);
    collider.arrow(0, 1), collider.arrow(2, 1), collider.arrow(1, 3);
    s.require("a chain a → m → b: open, and blocked by m", !pg::d_separated(chain, {0}, {2}, {}) && pg::d_separated(chain, {0}, {2}, {1}));
    s.require("a fork a ← m → b: open, and blocked by m", !pg::d_separated(fork, {0}, {2}, {}) && pg::d_separated(fork, {0}, {2}, {1}));
    s.require("a collider a → m ← b: blocked, opened by m and by a descendant of m",
              pg::d_separated(collider, {0}, {2}, {}) && !pg::d_separated(collider, {0}, {2}, {1}) && !pg::d_separated(collider, {0}, {2}, {3}));
  }
  {
    // M-bias: 0 = X, 1 = Y, 2 = U1, 3 = U2, 4 = M.
    pg::Dag g(5);
    g.arrow(2, 0), g.arrow(2, 4), g.arrow(3, 4), g.arrow(3, 1), g.arrow(0, 1);
    s.require("M-bias: the empty set satisfies the back-door criterion", pg::backdoor(g, 0, 1, {}));
    s.require("... {M} does not (it opens the collider), {M, U₁} does", !pg::backdoor(g, 0, 1, {4}) && pg::backdoor(g, 0, 1, {4, 2}));
    const auto sets = pg::adjustment_sets(g, 0, 1, {2, 3, 4});
    s.require("... the minimal sets among U₁, U₂, M: the empty set alone", sets.size() == 1 && sets[0].empty(), std::to_string(sets.size()) + " sets");
  }
  {
    // The front door: 0 = X, 1 = Y, 2 = M, 3 = U (unobserved).
    pg::Dag g(4);
    g.arrow(3, 0), g.arrow(3, 1), g.arrow(0, 2), g.arrow(2, 1);
    s.require("front door: {M} satisfies the front-door criterion", pg::frontdoor(g, 0, 1, {2}));
    s.require("... and no set of the observed (M) the back-door one", pg::adjustment_sets(g, 0, 1, {2}).empty());
    s.require("... with U observed, {U} is the back-door set", pg::adjustment_sets(g, 0, 1, {2, 3}).size() == 1 && pg::adjustment_sets(g, 0, 1, {2, 3})[0] == std::set<std::size_t>{3});
  }
  {
    calx::random::Stream st(20261108, 0x64616773u, 0, 0);  // "dags"
    int agree = 0, total = 0, separated = 0;
    for (int r = 0; r < 12; r++) {
      const std::size_t n = 7;
      pg::Dag g(n);
      std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.0));  // B[child][parent]
      for (std::size_t j = 1; j < n; j++)
        for (std::size_t i = 0; i < j; i++)
          if (st.uniform() < 0.35) {
            g.arrow(i, j);
            B[j][i] = (st.uniform() < 0.5 ? -1 : 1) * st.uniform(0.5, 1.5);
          }
      // X = B X + ε, ε ~ N(0, I): Σ = (I − B)⁻¹ (I − B)⁻ᵀ.
      calx::linalg::Matrix IB(n), I(n);
      for (std::size_t a = 0; a < n; a++) {
        I(a, a) = 1;
        for (std::size_t b = 0; b < n; b++) IB(a, b) = (a == b ? 1.0 : 0.0) - B[a][b];
      }
      const auto inv = calx::linalg::solve(IB, I);
      std::vector<std::vector<double>> S(n, std::vector<double>(n, 0.0));
      for (std::size_t a = 0; a < n; a++)
        for (std::size_t b = 0; b < n; b++)
          for (std::size_t k = 0; k < n; k++) S[a][b] += inv(a, k) * inv(b, k);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = i + 1; j < n; j++) {
          std::vector<std::size_t> rest;
          for (std::size_t k = 0; k < n; k++)
            if (k != i && k != j) rest.push_back(k);
          std::vector<std::set<std::size_t>> conds{{}};
          for (std::size_t a = 0; a < rest.size(); a++) {
            conds.push_back({rest[a]});
            for (std::size_t b = a + 1; b < rest.size(); b++) conds.push_back({rest[a], rest[b]});
          }
          for (const auto& zs : conds) {
            const bool d = pg::d_separated(g, {i}, {j}, zs);
            const bool zero = std::fabs(partial_correlation(S, i, j, zs)) < 1e-9;
            agree += d == zero, total++, separated += d;
          }
        }
    }
    s.require("random graphs: d-separation holds exactly where the linear Gaussian model's partial correlation vanishes", agree == total,
              std::to_string(agree) + " of " + std::to_string(total) + " (" + std::to_string(separated) + " separated)");
  }
  {
    pg::Dag cyc(3);
    cyc.arrow(0, 1), cyc.arrow(1, 2), cyc.arrow(2, 0);
    pg::Dag ok(3);
    ok.arrow(0, 1), ok.arrow(0, 2);
    s.require("acyclicity: a directed triangle is cyclic, a fork is not", !pg::acyclic(cyc) && pg::acyclic(ok));
  }
  return s.done();
}
