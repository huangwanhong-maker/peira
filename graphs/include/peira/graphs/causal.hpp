// Causal graphs: which conditional independences a directed acyclic graph implies, and which sets of variables identify an effect
// by adjustment (Pearl 1988; Pearl 2009, Causality, 2nd ed., ch. 3; Koller and Friedman 2009, Probabilistic Graphical Models, ch. 3).
//
// D-SEPARATION. A path between X and Y is blocked by Z when it holds a chain or fork a → m → b, a ← m → b with m ∈ Z, or a collider
// a → m ← b with neither m nor any descendant of m in Z. X and Y are d-separated by Z when every path between them is blocked; then
// X ⊥ Y | Z in every distribution that factorises over the graph, and in a distribution faithful to it only then (Pearl 2009,
// Thm 1.2.4; Spirtes, Glymour and Scheines 2000). Decided by the reachability of active trails ("Bayes ball"; Shachter 1998;
// Koller and Friedman 2009, Alg. 3.1), in O(V + E).
//
// THE BACK-DOOR CRITERION (Pearl 2009, Def. 3.3.1 and Thm 3.3.2). A set Z of observed variables satisfies it for (X, Y) when no
// element of Z descends from X, and Z blocks every path between X and Y that begins with an arrow into X: equivalently, X and Y are
// d-separated by Z in the graph with the arrows out of X removed. Then P(y | do(x)) = Σ_z P(y | x, z) P(z): adjustment identifies
// the effect.
//
// THE FRONT-DOOR CRITERION (Pearl 2009, Def. 3.3.3 and Thm 3.3.4). Z satisfies it for (X, Y) when (i) Z intercepts every directed
// path from X to Y, (ii) no back-door path from X to Z is unblocked (X ⊥ Z in the graph without X's outgoing arrows), and (iii) X
// blocks every back-door path from Z to Y. Then the effect is identified even when X and Y share an unobserved cause.
//
// ADJUSTMENT SETS. The minimal back-door sets among the observed non-descendants of X, by enumeration in order of size: for the
// small graphs a study draws (it is exponential in their number).

#pragma once

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace peira::graphs {

struct Dag {
  std::vector<std::vector<std::size_t>> parents, children;
  std::vector<std::string> names;
  explicit Dag(std::size_t n = 0) : parents(n), children(n), names(n) {}
  std::size_t size() const { return parents.size(); }
  void arrow(std::size_t from, std::size_t to) {
    if (from >= size() || to >= size() || from == to) throw std::invalid_argument("peira::graphs: an arrow between two distinct vertices");
    parents[to].push_back(from);
    children[from].push_back(to);
  }
};

// Whether the graph has no directed cycle (Kahn's algorithm).
inline bool acyclic(const Dag& g) {
  std::vector<std::size_t> in(g.size());
  for (std::size_t v = 0; v < g.size(); v++) in[v] = g.parents[v].size();
  std::vector<std::size_t> ready;
  for (std::size_t v = 0; v < g.size(); v++)
    if (in[v] == 0) ready.push_back(v);
  std::size_t seen = 0;
  while (!ready.empty()) {
    const std::size_t v = ready.back();
    ready.pop_back(), seen++;
    for (const std::size_t c : g.children[v])
      if (--in[c] == 0) ready.push_back(c);
  }
  return seen == g.size();
}

inline std::set<std::size_t> descendants(const Dag& g, std::size_t x) {
  std::set<std::size_t> out;
  std::vector<std::size_t> stack{x};
  while (!stack.empty()) {
    const std::size_t v = stack.back();
    stack.pop_back();
    for (const std::size_t c : g.children[v])
      if (out.insert(c).second) stack.push_back(c);
  }
  return out;
}

// Whether every vertex of xs is d-separated from every vertex of ys by zs.
inline bool d_separated(const Dag& g, const std::set<std::size_t>& xs, const std::set<std::size_t>& ys, const std::set<std::size_t>& zs) {
  for (const std::size_t v : xs)
    if (ys.count(v) || zs.count(v)) throw std::invalid_argument("peira::graphs: X, Y and Z disjoint");
  // The ancestors of Z (and Z): a collider is open when it is among them.
  std::set<std::size_t> anc(zs.begin(), zs.end());
  std::vector<std::size_t> stack(zs.begin(), zs.end());
  while (!stack.empty()) {
    const std::size_t v = stack.back();
    stack.pop_back();
    for (const std::size_t p : g.parents[v])
      if (anc.insert(p).second) stack.push_back(p);
  }
  // Reachability over (vertex, direction): up = reached from a child, down = reached from a parent (Koller and Friedman, Alg. 3.1).
  std::set<std::pair<std::size_t, bool>> visited;
  std::vector<std::pair<std::size_t, bool>> todo;
  for (const std::size_t x : xs) todo.push_back({x, true});
  while (!todo.empty()) {
    const auto [v, up] = todo.back();
    todo.pop_back();
    if (!visited.insert({v, up}).second) continue;
    if (!zs.count(v) && ys.count(v)) return false;
    if (up && !zs.count(v)) {
      for (const std::size_t p : g.parents[v]) todo.push_back({p, true});
      for (const std::size_t c : g.children[v]) todo.push_back({c, false});
    } else if (!up) {
      if (!zs.count(v))
        for (const std::size_t c : g.children[v]) todo.push_back({c, false});
      if (anc.count(v))
        for (const std::size_t p : g.parents[v]) todo.push_back({p, true});
    }
  }
  return true;
}

// The graph without the arrows out of x.
inline Dag without_outgoing(const Dag& g, std::size_t x) {
  Dag h(g.size());
  h.names = g.names;
  for (std::size_t v = 0; v < g.size(); v++)
    for (const std::size_t c : g.children[v])
      if (v != x) h.arrow(v, c);
  return h;
}

inline bool backdoor(const Dag& g, std::size_t x, std::size_t y, const std::set<std::size_t>& zs) {
  const auto desc = descendants(g, x);
  for (const std::size_t z : zs)
    if (desc.count(z) || z == x || z == y) return false;
  return d_separated(without_outgoing(g, x), {x}, {y}, zs);
}

inline bool frontdoor(const Dag& g, std::size_t x, std::size_t y, const std::set<std::size_t>& zs) {
  if (zs.count(x) || zs.count(y)) return false;
  // (i) no directed path from x to y avoiding Z.
  std::set<std::size_t> seen{x};
  std::vector<std::size_t> stack{x};
  while (!stack.empty()) {
    const std::size_t v = stack.back();
    stack.pop_back();
    for (const std::size_t c : g.children[v]) {
      if (zs.count(c)) continue;
      if (c == y) return false;
      if (seen.insert(c).second) stack.push_back(c);
    }
  }
  const Dag gx = without_outgoing(g, x);
  for (const std::size_t z : zs) {
    if (!d_separated(gx, {x}, {z}, {})) return false;                     // (ii)
    if (!d_separated(without_outgoing(g, z), {z}, {y}, {x})) return false;  // (iii)
  }
  return true;
}

// The minimal back-door sets for (x, y) among the vertices `observed` that do not descend from x.
inline std::vector<std::set<std::size_t>> adjustment_sets(const Dag& g, std::size_t x, std::size_t y, const std::vector<std::size_t>& observed) {
  const auto desc = descendants(g, x);
  std::vector<std::size_t> pool;
  for (const std::size_t v : observed)
    if (v != x && v != y && !desc.count(v)) pool.push_back(v);
  if (pool.size() > 20) throw std::invalid_argument("peira::graphs: adjustment sets are enumerated among 20 candidates at most");
  std::vector<std::set<std::size_t>> found;
  for (std::size_t size = 0; size <= pool.size(); size++) {
    std::vector<char> pick(pool.size(), 0);
    std::fill(pick.end() - static_cast<std::ptrdiff_t>(size), pick.end(), 1);
    do {
      std::set<std::size_t> zs;
      for (std::size_t k = 0; k < pool.size(); k++)
        if (pick[k]) zs.insert(pool[k]);
      const bool contains_found = std::any_of(found.begin(), found.end(), [&](const std::set<std::size_t>& f) { return std::includes(zs.begin(), zs.end(), f.begin(), f.end()); });
      if (!contains_found && backdoor(g, x, y, zs)) found.push_back(zs);
    } while (std::next_permutation(pick.begin(), pick.end()));
  }
  return found;
}

}  // namespace peira::graphs
