# 0006. Interference: Treatments That Reach Beyond Their Units

- **Status:** Accepted, provisionally (2026-11-09).
- **Decided by:** the infrastructure's ADR 0087, the urban priorities' stage 3, at the user's word ("Good! Begin.") on the proposal that the evaluation of street changes be built as general layers ("does this module / algorithm can be reused? We don't want it just for a city"); the module's content, Claude.

## 1. Context

Every estimator in `quasi` and `randomisation` assumes that a unit's outcome depends on its own treatment alone. Treatments in space break it:
- a closed street sends walkers to the next;
- a cycle lane draws riders from parallel streets;
- a low-emission zone cleans the air beyond its edge;
- a vaccinated child protects a classmate.

The spillover is an effect of its own, and it contaminates the comparison units near the treatment. The methods that restore identification name no object in nature. They are designs of exposure and what they identify, so they are peira's.

## 2. Decision

A module `interference` (§9 of the theory note):
- `exposure.hpp`: an exposure mapping by rings of distance (Aronow and Samii 2017) from treatments that begin at different times, with each unit's ring by period, its onsets and its entry ring; the distances are the caller's, in any metric;
- `rings.hpp`:
  - the two-way regression on ring indicators (Butts 2021);
  - the staggered estimator by ring, against units not yet reached;
  - the clean donors for a synthetic control;
- `spillover_test.hpp`: Athey, Eckles and Imbens's exact test of no spillovers in a completely randomised design.

`quasi`'s staggered estimator gains an optional `reached`, each unit's first period touched by any treatment. A unit is a comparison only while not yet reached. Without it, the estimator is unchanged.

## 3. Alternatives Considered

- **In the domains.** Rejected: urban science, epidemiology and ecology would each write the same mapping.
- **A general interference graph** (Hudgens and Halloran's partial interference, network exposure with arbitrary mappings). Deferred: rings of distance cover the spatial cases asked for, and the mapping is one function, so others can join beside it.
- **Modelling the spillover inside the synthetic control** (Cao and Dowd 2019). Deferred: leaving the reached donors out is what the mapping justifies, and is simpler to state.

## 4. Verification

`check_peira_interference` (14): exposure by hand; both effects recovered to rounding without noise, and within four standard errors with it; the bias of ignoring spillovers to $10^{-13}$; ring-staggered effects exact with zero placebos; the spilled-over donor's bias as $4 - 3w$ exactly; the exactness of the test by enumeration, and its power. `quasi`'s checks pass unchanged.

## 5. Consequences

- The infrastructure's urban science evaluates street changes on any city's sensors through it. Other sciences can reach it too: epidemiology (planned) and ecology.
- peira's checks rise from 146 to 160.
