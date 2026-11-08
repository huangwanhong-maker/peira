# Survival

Status: **draft** (2026-11-08). Module `peira::survival`.

## 1. Concepts

Times to an event, seen through censoring, built from risk sets and the events at each time. The theory, with proofs, is in [docs/theory.md](../docs/theory.md).

## 2. Content

- `kaplan_meier` with Greenwood's variance and log(−log) intervals, its median; `survival_at`.
- `nelson_aalen`.
- `log_rank` for k groups.
- `cox`: proportional hazards with Breslow's or Efron's ties, Newton–Raphson, the baseline cumulative hazard.
- `aalen_johansen`: cumulative incidence with competing causes.

## 3. Checks

statsmodels 0.15.0 and lifelines 0.30.3 on Rossi's recidivism study and Walton's lifetimes (`tests/fixtures/make_survival_fixtures.py`): every estimate, standard error and interval to about 10⁻¹⁴; exact identities (Breslow = Efron without ties; the score test is the log-rank; incidences sum to 1 − S).

## 4. Dependencies

tyche `probability`, calx `linalg`.

## 5. References

Kaplan and Meier (1958); Greenwood (1926); Nelson (1972), Aalen (1978); Mantel (1966); Cox (1972); Breslow (1974); Efron (1977); Aalen and Johansen (1978); Kalbfleisch and Prentice (2002).
