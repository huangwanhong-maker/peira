# 0007. Stated-Choice Designs

- **Status:** Accepted, provisionally (2026-11-09).
- **Decided by:** the infrastructure's ADR 0088 (the urban priorities' stage 4, demand), at the user's word ("continue the next feature").

## 1. Decision

`design` gains `choice.hpp`, the design of stated-choice experiments:
- attributes categorical (effects coding) or numeric;
- the information of a design under the logit at priors, and its D-error, utility-neutral, at priors or Bayesian over draws (Huber and Zwerina 1996; Sándor and Wedel 2001);
- coordinate exchange from random starts, never keeping a set with two equal profiles (Meyer and Nachtsheim 1995; Kessels, Jones, Goos and Vandebroek 2009);
- level counts.

## 2. Rationale

A stated-choice survey is an experiment whose design the experimenter chooses: which alternatives to show together decides what the answers identify. It is a design of the mechanism of observation, so it is peira's. Estimation is nomos's conditional logit.

## 3. Verification

`check_peira_choice_designs` (8):
- **The bound, exactly.** For binary attributes in pairs, D₀ ≥ 1/S, attained by the foldover of an orthogonal array: D₀ = 1/4 with information S I.
- **The search reaches it.** Coordinate exchange finds D₀ = 1/8 for 8 pairs.
- **A pair's information** is ddᵀ/4.
- **Priors.** A design made at priors is at least as good there as the utility-neutral one.

peira's checks rise from 160 to 168.
