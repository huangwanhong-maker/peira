"""Reference values for peira/survival, made by running statsmodels and lifelines (never written from memory).

Data, as lifelines 0.30 ships them (lifelines.datasets):
- rossi: Rossi, Berk and Lenihan (1980), Money, Work and Crime: 432 men released from Maryland state prisons, followed for a
  year; week of first arrest (arrest = 1) or censoring at week 52, with financial aid (fin), age, race, work experience
  (wexp), marital status (mar), parole (paro) and prior convictions (prio). The standard example of Fox (2002), An R and
  S-PLUS Companion to Applied Regression, appendix on Cox regression.
- waltons: lifetimes of 163 flies (T, E) in two groups, miR-137 and control, from lifelines' documentation.

Values written:
- the Kaplan-Meier estimate at each event time with Greenwood's standard error (statsmodels SurvfuncRight) and the
  log(-log) interval (lifelines KaplanMeierFitter, its default);
- the Nelson-Aalen cumulative hazard without tie smoothing (lifelines NelsonAalenFitter(nelson_aalen_smoothing=False));
- the log-rank test (statsmodels survdiff; lifelines logrank_test as a cross-check);
- Cox's proportional hazards on rossi with Breslow's and Efron's ties (statsmodels PHReg), and Efron's by lifelines
  (CoxPHFitter) as a cross-check.

Run with an environment holding statsmodels and lifelines:  python make_survival_fixtures.py
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path

import lifelines
import numpy as np
import statsmodels
from lifelines import CoxPHFitter, KaplanMeierFitter, NelsonAalenFitter
from lifelines.datasets import load_rossi, load_waltons
from lifelines.statistics import logrank_test
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

HERE = Path(__file__).resolve().parent
COVARIATES = ["fin", "age", "race", "wexp", "mar", "paro", "prio"]

rossi = load_rossi()
waltons = load_waltons()
waltons["g"] = (waltons["group"] == "miR-137").astype(int)

rossi[["week", "arrest", *COVARIATES]].to_csv(HERE / "rossi.csv", index=False)
waltons[["T", "E", "g"]].to_csv(HERE / "waltons.csv", index=False)

out = {
    "made": datetime.date.today().isoformat(),
    "by": f"statsmodels {statsmodels.__version__}, lifelines {lifelines.__version__}, numpy {np.__version__}",
}

# Kaplan-Meier on rossi.
sf = SurvfuncRight(rossi["week"].to_numpy(), rossi["arrest"].to_numpy())
km = KaplanMeierFitter().fit(rossi["week"], rossi["arrest"])
ci = km.confidence_interval_
times = [float(t) for t in sf.surv_times]
lower = [float(ci.loc[t].iloc[0]) for t in times]
upper = [float(ci.loc[t].iloc[1]) for t in times]
out["km_rossi"] = {"time": times, "survival": [float(v) for v in sf.surv_prob], "se": [float(v) for v in sf.surv_prob_se],
                   "lower": lower, "upper": upper}

# Kaplan-Meier on waltons, whose median is reached.
km_w = KaplanMeierFitter().fit(waltons["T"], waltons["E"])
out["km_waltons_median"] = float(km_w.median_survival_time_)

# Nelson-Aalen on rossi, without tie smoothing.
na = NelsonAalenFitter(nelson_aalen_smoothing=False).fit(rossi["week"], rossi["arrest"])
out["na_rossi"] = {"time": times, "hazard": [float(na.cumulative_hazard_.loc[t].iloc[0]) for t in times]}

# Log-rank.
for name, frame, t, e, g in (("logrank_waltons", waltons, "T", "E", "g"), ("logrank_rossi_fin", rossi, "week", "arrest", "fin")):
    chi2, p = survdiff(frame[t].to_numpy(), frame[e].to_numpy(), frame[g].to_numpy())
    a, b = frame[frame[g] == 1], frame[frame[g] == 0]
    lr = logrank_test(a[t], b[t], a[e], b[e])
    out[name] = {"statistic": float(chi2), "p_value": float(p), "lifelines_statistic": float(lr.test_statistic), "lifelines_p_value": float(lr.p_value)}

# Cox on rossi.
X = rossi[COVARIATES].to_numpy(dtype=float)
for ties in ("breslow", "efron"):
    fit = PHReg(rossi["week"].to_numpy(dtype=float), X, status=rossi["arrest"].to_numpy(), ties=ties).fit()
    out[f"cox_rossi_{ties}"] = {"beta": [float(v) for v in fit.params], "se": [float(v) for v in fit.bse], "log_partial_likelihood": float(fit.llf)}
cph = CoxPHFitter().fit(rossi[["week", "arrest", *COVARIATES]], duration_col="week", event_col="arrest")
out["cox_rossi_efron_lifelines"] = {"beta": [float(cph.params_[c]) for c in COVARIATES], "se": [float(cph.standard_errors_[c]) for c in COVARIATES],
                                    "log_partial_likelihood": float(cph.log_likelihood_)}

(HERE / "values.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


# The same values as lines "name v1 v2 ...", which the C++ checks read without a JSON parser.
def flat(prefix, value, lines):
    if isinstance(value, dict):
        for k, v in value.items():
            flat(f"{prefix}.{k}" if prefix else k, v, lines)
    elif isinstance(value, list):
        lines.append(prefix + " " + " ".join(repr(float(v)) for v in value))
    elif isinstance(value, (int, float)):
        lines.append(f"{prefix} {float(value)!r}")


lines = [f"# made {out['made']} by {out['by']} (make_survival_fixtures.py)"]
for key, value in out.items():
    if key not in ("made", "by"):
        flat(key, value, lines)
(HERE / "values.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("written:", ", ".join(k for k in out if k not in ("made", "by")))
