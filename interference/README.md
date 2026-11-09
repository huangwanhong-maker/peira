# Interference

Status: **draft** (2026-11-09). Module `peira::interference`, [decision 0006](../development/decisions/0006-interference.md).

## 1. Concepts

A treatment given to one unit may reach others: a closed street sends its walkers to the next, a vaccinated child protects a classmate, a low-emission zone cleans the air beyond its edge, a fertilised plot feeds the plot beside it. The stable unit treatment value assumption, on which every estimator in `quasi` and `randomisation` rests, then fails. Interference does two things:
- **it is an effect of its own,** the spillover;
- **it contaminates the comparison,** since the untreated units near a treatment are no longer untreated.

An exposure mapping (Aronow and Samii 2017) restores identification by an assumption that can be stated. A unit's outcome depends on the whole assignment only through its ring of distance from the nearest treatment begun, and nothing reaches beyond the outermost ring. Distances are the caller's, in whatever metric carries the effect: metres along a network, kilometres over the ground, steps in a graph. The theory is in [docs/theory.md](../docs/theory.md), §9.

## 2. Content

- `exposure.hpp`:
  - `Contact`, `expose`: each unit's ring in each period, its onsets by ring and its entry ring, from treatments that begin at different times;
  - `rings_of`: the ring of each observation of a long panel.
- `rings.hpp`:
  - `ring_two_way`: the two-way regression on ring indicators, giving the direct effect and the spillover in each ring (Butts 2021), with errors clustered by unit;
  - `ring_staggered`: the group-time effects of entering each ring, against units not yet reached (Callaway and Sant'Anna 2021 by ring), through quasi's staggered estimator with `reached`;
  - `clean_donors`: the units never reached, the donors a synthetic control may use.
- `spillover_test.hpp`:
  - `no_spillover`: Athey, Eckles and Imbens's exact randomisation test of no spillovers in a completely randomised design, conditional on the focal units' treatments.
- `interference.hpp`: all three.

## 3. Checks

`check_peira_interference` (14):
- **Exposure.** A schedule of two treatments by hand: rings, onsets, entry rings.
- **Rings in a two-way regression.**
  - Without noise, the direct effect and the spillover recovered to rounding.
  - The regression that ignores spillovers is off by exactly $n_s/(n_s + n_u)\,\tau_1$.
  - With noise, each effect within four standard errors.
- **Rings with timing that varies.**
  - Without noise, both effects recovered exactly, and the placebo periods are zero.
  - The plain staggered estimator, counting spilled-over units as comparisons, is biased.
- **Synthetic control.**
  - A spilled-over twin among the donors takes weight $w$ and biases the gap to $4 - 3w$ exactly.
  - The clean donors give the effect.
- **The test of no spillovers.**
  - Exact under the null in every conditioning class, by enumeration.
  - It finds a spillover of 3 on 24 units, with $p = 0.002$.

## 4. Dependencies

peira `quasi` (two-way regressions, staggered adoption, synthetic control) and `design` (enumeration and draws of complete randomisation); calx `linalg` and `random`.

## 5. References

- Aronow, P. M. and Samii, C. (2017), Estimating average causal effects under general interference, with application to a social network experiment, *Ann. Appl. Stat.* 11, 1912–1947.
- Athey, S., Eckles, D. and Imbens, G. W. (2018), Exact p-values for network interference, *J. Am. Stat. Assoc.* 113, 230–240.
- Butts, K. (2021), Difference-in-differences estimation with spatial spillovers, arXiv:2105.03737.
- Callaway, B. and Sant'Anna, P. H. C. (2021), Difference-in-differences with multiple time periods, *J. Econometrics* 225, 200–230.
- Cao, J. and Dowd, C. (2019), Estimation and inference for synthetic control methods with spillover effects, arXiv:1902.07343.
- Clarke, D. (2017), Estimating difference-in-differences in the presence of spillovers, MPRA paper 81604.
- Cox, D. R. (1958), *Planning of Experiments*, Wiley.
- Hudgens, M. G. and Halloran, M. E. (2008), Toward causal inference with interference, *J. Am. Stat. Assoc.* 103, 832–842.
- Rubin, D. B. (1980), Comment on "Randomization analysis of experimental data: the Fisher randomization test" by D. Basu, *J. Am. Stat. Assoc.* 75, 591–593.
