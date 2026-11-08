# CLAUDE.md

This file guides Claude Code in the peira repository.

peira is one of the foundation libraries of [Noetarium, the Science and Humanity Infrastructure](https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure), which includes it as the git submodule `foundations/peira`. It is a project in its own right and must make sense without the infrastructure.

## 1. Purpose and Scope

peira implements the designs and mechanisms by which units come to be treated and observed, and what is identified through them:

- potential outcomes and estimands; randomised designs, their power and randomisation inference;
- weighting by the probability of being seen; designs of samples and their estimates;
- survival and censoring; cohorts;
- missing data and multiple imputation;
- quasi-experiments; adjustment sets on causal graphs;
- meta-analysis, sequential testing, sensitivity and partial identification, pre-registered plans.

It builds on four sibling libraries:

- **calx**, for random streams and linear algebra;
- **morphe**, for graphs (causal graphs, later);
- **physis**, as the siblings' order requires;
- **tyche**, for distributions, regression, generalised linear models and permutation p-values.

Out of scope:

| Topic | Where it belongs |
|---|---|
| Inference within a given model: likelihood, estimators, regression, panels with fixed effects, multiple testing, time series | tyche |
| Measures and mechanisms of social life: inequality, index numbers, games, matching markets, voting, agents | nomos |
| A particular population, disease, trial or city | the domains of the infrastructure |
| Drawing survival curves, forest plots, Lexis diagrams | the infrastructure's `visualization/` |

The boundary is [decision 0002](development/decisions/0002-scope-boundaries-and-dependencies.md): a method whose object is the mechanism of assignment, sampling, censoring or response, and what it identifies, is peira's.

## 2. Invariants

1. **Dependencies.** peira depends only on calx, morphe, physis and tyche among project libraries, and never on nomos or on the infrastructure's domains, edges or applications.
2. **No particular objects.** Nothing in peira names an object in nature.
3. **Pure.** No I/O, no global state. Every draw comes from a calx stream given by the caller, or seeded explicitly.
4. **Truth from the design.** Wherever a design can be enumerated, a check enumerates it and holds its expectations exactly. Otherwise reference values come from running a reference library (a script in `tests/fixtures/`, never from memory), or a simulation states its Monte Carlo error.
5. **Identification stated.** Each estimator's header says what it identifies and under which assumptions, and which of them cannot be tested.
6. **Reuse.** What a sibling holds is used through its public headers, never written again here.
7. **Cited sources.** Every implemented method cites its source; equation or section numbers are given only where certain.

## 3. Design Philosophy

peira shares its governing principle with its sibling repositories: **modular, extensible, revisable**.

- **Status vocabulary:**
  - `planned`: named, no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, with a pointer to its replacement.
- **One structure.** Experiments, surveys, censoring and nonresponse are full data coarsened by a mechanism ([docs/theory.md](docs/theory.md)). A method that applies to one of them is written so that it applies to all where it can: inverse-probability weighting is one header for samples and arms.

## 4. Repository Layout

```
peira/
├── CMakeLists.txt   the project: one header-only target peira::<module> per module (peira_module); the siblings added once
├── design/ · randomisation/ · weighting/ · surveys/ · survival/ · missing/
├── docs/            the theory note
└── development/     decisions · ideas · release notes
```

Each module directory holds:
- `README.md`: concepts first, then the content, then the checks;
- `include/peira/<module>/`;
- `tests/check_<module>.cpp`, a check program on calx's `verification`, registered with CTest by `peira_module` in the top-level `CMakeLists.txt`, with reference values in `tests/fixtures/` (`PEIRA_FIXTURES`).

The modules and their order are in [decision 0003](development/decisions/0003-the-first-modules.md).

## 5. Working Procedures

### 5.1 Decisions and Ideas

- Decisions are numbered records in `development/decisions/`. Records are never rewritten; a later record supersedes an earlier one.
- Ideas are files named `development/ideas/YYYY-MM-DD-slug.md`, with status `seed`, `promoted` or `abandoned`.
- **Release notes.** At each commit, add to `development/release-notes/YYYY-MM-DD.md` (one file per day, one short section per commit) a concise summary in plain prose of what is new, what one can now do and what works better. It is not a list of changes (the git log keeps those), and it is not written between commits. This convention is shared with the sibling repositories.

### 5.2 Change Discipline

- Before restructuring or changing an invariant (§2), propose the change and wait for the user's agreement.
- Keep each change to one concept.

### 5.3 Shared Conventions

calx, morphe, physis, tyche, peira, nomos, eikon, heuresis and the infrastructure keep these conventions identical:

- writing style;
- status vocabulary;
- decision, idea and release-note conventions.

Change them in every repository together.

### 5.4 Working Inside the Infrastructure

When you edit peira inside the infrastructure's `foundations/peira`:

1. commit in peira first;
2. then commit the updated pointer in the infrastructure.

When a change needs a sibling to change as well, commit them in dependency order: calx and morphe first, then physis, then tyche, then peira.

## 6. Conventions

- **Arms.** z = 1 is treated, z = 0 is control; effects are treated minus control.
- **Weights.** A weight is the inverse of a probability of being seen, w = 1/π.
- **Censoring.** event = 1 is an event seen; 0 is censored at the time given.
- **Writing style.**
  - Headings are declarative, numbered and academic, never questions.
  - Concepts come before structure.
  - Label asides with *Rationale.*, *Remark.* or *Caveat.*
  - British spelling; English throughout.

## 7. Technology

- **Language: C++.** This follows ADR 0003 of the infrastructure.
- **Provisional:** C++20, header-only targets `peira::<module>`; CMake 3.20 or later. Checks are programs on calx's `verification`, registered with CTest. On the development machine the compiler is GCC in WSL.
- **Dependency resolution.** Inside the infrastructure, peira uses the targets of `foundations/calx`, `foundations/morphe`, `foundations/physis` and `foundations/tyche`, adding each that no larger project has added. Never nest submodules.
