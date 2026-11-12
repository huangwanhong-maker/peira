// Imperfect detection (peira's decision 0008; docs/theory.md, §10): full data (whether a site is occupied, how many individuals a
// site or a population holds, where individuals lie) seen through a mechanism of detection that misses some of them, and the
// designs that identify it. Site occupancy from repeated visits (occupancy.hpp), abundance from repeated counts (n_mixture.hpp),
// closed populations from captures and recaptures (capture_recapture.hpp), density from distances (distance.hpp), on the
// likelihood of latent counts seen through binomial detection (latent.hpp); surveys drawn from known truths (simulate.hpp).

#pragma once

#include "capture_recapture.hpp"
#include "common.hpp"
#include "distance.hpp"
#include "latent.hpp"
#include "n_mixture.hpp"
#include "occupancy.hpp"
#include "simulate.hpp"
