// Surveys drawn from known full data through a known mechanism of detection: for planning a design (how many sites, visits,
// occasions or lines a precision needs) and for the checks. Every draw comes from the caller's calx stream.
//
// - simulate_visits: at each site a latent count from the mixing law (presence with probability ψ_i; a Poisson count of mean
//   λ_i; a negative binomial of mean λ_i and size r, drawn as a Poisson whose mean is Gamma(r, λ_i/r)), then on each visit
//   y_ij ~ Binomial(N_i, p_ij); a visit whose p_ij is NaN is not made (−1).
// - simulate_captures: N individuals, each caught on occasion j with probability p_j independently; the histories of those
//   caught at least once.
// - simulate_distances: individuals at density D. For lines, a Poisson number with mean 2wDl_k within w of line k, each at a
//   perpendicular distance uniform on [0, w]; for points, on each visit a Poisson number with mean πw²D within w of the point,
//   at radius w√U (uniform in the disc), the visits independent (individuals redistributed between them). Each is detected with
//   probability g(distance).

#pragma once

#include <calx/random/stream.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <span>
#include <vector>

#include "capture_recapture.hpp"
#include "common.hpp"
#include "distance.hpp"
#include "latent.hpp"

namespace peira::detection {

struct SimulatedVisits {
  SiteVisits y;
  std::vector<int> latent;  // N_i (z_i for presence): the full data
};

inline SimulatedVisits simulate_visits(calx::random::Stream& rng, Mixing mixing, std::span<const double> mean, std::span<const double> p, std::size_t visits,
                                       double size = 0) {
  const std::size_t S = mean.size();
  detail::require(p.size() == S * visits, "one detection probability per site and visit");
  detail::require(mixing != Mixing::negative_binomial || size > 0, "a size r > 0 for the negative binomial");
  SimulatedVisits out;
  out.y = SiteVisits(S, visits);
  out.latent.resize(S);
  for (std::size_t i = 0; i < S; i++) {
    std::int64_t N = 0;
    switch (mixing) {
      case Mixing::bernoulli: N = rng.bernoulli(mean[i]) ? 1 : 0; break;
      case Mixing::poisson: N = rng.poisson(mean[i]); break;
      default: N = rng.poisson(rng.gamma(size, mean[i] / size)); break;
    }
    out.latent[i] = static_cast<int>(N);
    for (std::size_t j = 0; j < visits; j++) {
      const double pij = p[i * visits + j];
      if (std::isnan(pij)) continue;
      out.y(i, j) = static_cast<int>(rng.binomial(N, pij));
    }
  }
  return out;
}

inline CaptureHistories simulate_captures(calx::random::Stream& rng, std::size_t N, std::span<const double> p) {
  const std::size_t t = p.size();
  std::vector<int> rows;
  std::size_t caught = 0;
  std::vector<int> h(t);
  for (std::size_t i = 0; i < N; i++) {
    bool any = false;
    for (std::size_t j = 0; j < t; j++) any |= (h[j] = rng.bernoulli(p[j]) ? 1 : 0);
    if (!any) continue;
    rows.insert(rows.end(), h.begin(), h.end());
    caught++;
  }
  CaptureHistories out(caught, t);
  out.y = std::move(rows);
  return out;
}

template <class G>
inline DistanceSurvey simulate_distances(calx::random::Stream& rng, Transect transect, double density, std::span<const double> effort, double w, G&& g) {
  detail::require(density >= 0 && w > 0 && std::isfinite(w), "a density ≥ 0 and a finite w > 0");
  DistanceSurvey s;
  s.transect = transect;
  s.truncation = w;
  s.effort.assign(effort.begin(), effort.end());
  for (std::size_t k = 0; k < effort.size(); k++) {
    if (transect == Transect::line) {
      const std::int64_t m = rng.poisson(2 * w * density * effort[k]);
      for (std::int64_t a = 0; a < m; a++) {
        const double x = w * rng.uniform();
        if (rng.bernoulli(g(x))) s.distance.push_back(x), s.on.push_back(k);
      }
    } else {
      const auto visits = static_cast<std::int64_t>(std::llround(effort[k]));
      for (std::int64_t v = 0; v < visits; v++) {
        const std::int64_t m = rng.poisson(std::numbers::pi * w * w * density);
        for (std::int64_t a = 0; a < m; a++) {
          const double r = w * std::sqrt(rng.uniform());
          if (rng.bernoulli(g(r))) s.distance.push_back(r), s.on.push_back(k);
        }
      }
    }
  }
  return s;
}

}  // namespace peira::detection
