# The Theory of peira: Full Data, Coarsening Mechanisms and Identification

*Status: draft (2026-11-08). The results the code relies on, with their hypotheses and proofs or citations. Each section names the module that implements it.*

## 1. The Coarsened-Data Model

**Definition 1.1 (full data and coarsening).** Let $U = \{1, \dots, N\}$ be a finite population or a sample of independent units. Unit $i$ carries *full data* $F_i$, which are never wholly seen, and a *coarsening variable* $R_i$. What is observed is $O_i = \Phi(F_i, R_i)$, a function that reveals part of $F_i$ according to $R_i$. The law of $R$ given $F$ is the *mechanism*.

| Setting | Full data $F_i$ | Coarsening $R_i$ | Observed $O_i$ |
|---|---|---|---|
| Experiment | $(Y_i(0), Y_i(1), X_i)$ | the arm $Z_i \in \{0, 1\}$ | $(Z_i, Y_i(Z_i), X_i)$ |
| Survey | $y_i$ for every $i \in U$ | inclusion $I_i = [i \in s]$ | $y_i$ for $i \in s$ |
| Censoring | event time $T_i$, covariates | censoring time $C_i$ | $(\min(T_i, C_i), [T_i \le C_i])$ |
| Nonresponse | $(Y_i, X_i)$ | response $R_i$ | $(X_i, R_i, R_i Y_i)$ |

**Definition 1.2 (identification).** A functional $\theta(P_F)$ of the full data's law is *identified* under a set of assumptions $\mathcal{A}$ on the mechanism when any two full-data laws compatible with $\mathcal{A}$ that induce the same law of $O$ have the same $\theta$.

*Remark.* Every theorem below has this shape. An assumption on the mechanism (it is known, or it ignores the hidden part) makes the observed law determine the estimand. Estimation then follows the identifying formula.

**Definition 1.3 (coarsening at random).** The mechanism is *coarsening at random* (CAR) when $P(R = r \mid F)$ depends on $F$ only through what $r$ reveals of it (Heitjan and Rubin 1991). Randomised assignment, probability sampling, independent censoring and MAR nonresponse are its cases.

## 2. Weighting by the Probability of Being Seen

*Module `weighting`, and through it `surveys`.*

**Theorem 2.1 (Horvitz and Thompson 1952).** Let a design on $U$ give each unit an inclusion probability $\pi_i = P(I_i = 1) > 0$. Then $\hat t = \sum_{i \in s} y_i/\pi_i$ satisfies $E[\hat t] = \sum_{i \in U} y_i$.

*Proof.* $\hat t = \sum_{i \in U} I_i y_i / \pi_i$, and $E[I_i] = \pi_i$. ∎

**Theorem 2.2 (the variance and its unbiased estimators).** With joint inclusion probabilities $\pi_{ij} = P(I_i I_j = 1)$, $\pi_{ii} = \pi_i$,
$$\operatorname{Var}(\hat t) = \sum_{i, j \in U} (\pi_{ij} - \pi_i \pi_j) \frac{y_i}{\pi_i} \frac{y_j}{\pi_j},$$
and if every $\pi_{ij} > 0$ then $\hat V_{HT} = \sum_{i, j \in s} \frac{\pi_{ij} - \pi_i\pi_j}{\pi_{ij}} \frac{y_i}{\pi_i}\frac{y_j}{\pi_j}$ is unbiased for it. For designs of fixed size, so is the Sen–Yates–Grundy form $\hat V_{SYG} = \tfrac12 \sum_{i \ne j \in s} \frac{\pi_i\pi_j - \pi_{ij}}{\pi_{ij}} (y_i/\pi_i - y_j/\pi_j)^2$.

*Proof.* The variance is $\sum_{i,j} \operatorname{Cov}(I_i, I_j) \frac{y_i y_j}{\pi_i \pi_j}$, with $\operatorname{Cov}(I_i, I_j) = \pi_{ij} - \pi_i\pi_j$. Each term of $\hat V_{HT}$ appears with probability $\pi_{ij}$ and is divided by it. For a fixed size $n$, $\sum_j (\pi_{ij} - \pi_i\pi_j) = 0$ for every $i$, which turns one form into the other (Särndal, Swensson and Wretman 1992, Result 2.8.1). ∎

*Checked* by enumerating every sample of an unequal-probability design of size 2 on five units: the expectations agree to $10^{-13}$.

**Corollary 2.3 (the effect as two totals).** Under a Bernoulli design with $e_i = P(Z_i = 1) \in (0, 1)$,
$$\hat\tau_{IPW} = \frac1N \sum_i \Big( \frac{Z_i Y_i}{e_i} - \frac{(1 - Z_i) Y_i}{1 - e_i} \Big)$$
is unbiased for $\tau = \frac1N \sum_i (Y_i(1) - Y_i(0))$. Each arm is a Horvitz–Thompson total, with the arm as the sample and $e_i$, $1 - e_i$ as the inclusion probabilities. *Checked* over all $2^6$ assignments, exactly.

**Theorem 2.4 (double robustness; Robins, Rotnitzky and Zhao 1994; Scharfstein, Rotnitzky and Robins 1999).** Under unconfoundedness $Z \perp (Y(0), Y(1)) \mid X$ and overlap $0 < e(X) < 1$, let
$$\psi = m_1(X) - m_0(X) + \frac{Z (Y - m_1(X))}{e(X)} - \frac{(1 - Z)(Y - m_0(X))}{1 - e(X)}.$$
Then $E[\psi] = \tau$ if either $e$ is the true propensity or $m_z(x) = E[Y \mid X = x, Z = z]$ for $z = 0, 1$.

*Proof.* Write $E[\psi] = E[m_1 - m_0] + E\big[\frac{Z(Y(1) - m_1)}{e}\big] - E\big[\frac{(1-Z)(Y(0) - m_0)}{1-e}\big]$.
- If $e$ is true, condition on $X$. By unconfoundedness $E[Z (Y(1) - m_1) / e \mid X] = E[Y(1) - m_1 \mid X]$, so the first two terms give $E[Y(1)]$, and likewise for the arm of controls.
- If the $m_z$ are true, $E[Z(Y(1) - m_1) \mid X] = e^\star(X)\, E[Y(1) - m_1(X) \mid X] = 0$ for whatever $e$ is used, so $E[\psi] = E[m_1 - m_0] = \tau$. ∎

$\psi$ is the efficient influence function of $\tau$ in the nonparametric model (Hahn 1998). Its sample variance over $n$ gives the standard error. *Checked* by simulation: unbiased when either nuisance is right, biased (30 Monte Carlo standard errors) when both are wrong.

## 3. Randomised Experiments

*Module `randomisation`.*

**Theorem 3.1 (Neyman 1923).** Fix potential outcomes on $n$ units and assign $n_1$ to treatment by complete randomisation. Then the difference in means $\hat\tau$ is unbiased for $\tau$, and
$$\operatorname{Var}(\hat\tau) = \frac{S_1^2}{n_1} + \frac{S_0^2}{n_0} - \frac{S_\tau^2}{n},$$
with $S_z^2$ the population variances of $Y(z)$ and $S_\tau^2$ that of the unit effects. The estimator $\hat V = s_1^2/n_1 + s_0^2/n_0$ has $E[\hat V] = S_1^2/n_1 + S_0^2/n_0$. It is conservative, and exact when the effect is constant.

*Proof.* Imbens and Rubin (2015), Theorem 6.2. The sample variance of each arm under simple random sampling without replacement is unbiased for its population variance. ∎ *Checked* over all 70 assignments of eight units, exactly, for a heterogeneous effect (strictly conservative) and a constant one (equal).

**Theorem 3.2 (Fisher's exact test; Fisher 1935).** Under the sharp null $Y_i(1) = Y_i(0) + \tau_0$ for every $i$, every potential outcome is known from the observed data. Hence the law of any statistic $T(Z)$ over the design is known, and the test rejecting when $p(Z) = P_{Z'}(T(Z') \text{ at least as extreme as } T(Z)) \le \alpha$ has level at most $\alpha$.

*Proof.* $p(Z)$ is the upper tail of $T$'s known law evaluated at $T(Z)$, and $P(p(Z) \le \alpha) \le \alpha$ for the tail function of any discrete law. ∎ *Checked* over all 252 assignments (12 of 252 rejected at 0.05), and on Fisher's lady tasting tea ($p = 1/70$).

## 4. Designs of Samples

*Module `surveys`.*

**Theorem 4.1 (the ultimate-cluster variance).** Let strata $h$ draw $n_h \ge 2$ primary units independently. Write $z_{ha}$ for the weighted total of a linearised variable over primary unit $a$. If the primary units are drawn with replacement (Hansen–Hurwitz), or if they are the units themselves drawn by simple random sampling without replacement at fraction $f_h$, then
$$\hat V = \sum_h (1 - f_h) \frac{n_h}{n_h - 1} \sum_{a \in h} (z_{ha} - \bar z_h)^2$$
is unbiased for the variance of the estimated total. Here $f_h = 0$ with replacement. Otherwise it is the standard approximation (Wolter 2007, §1.4; Särndal et al. 1992, Result 4.5.1).

*Checked* exactly in both cases by enumerating every sample. A ratio $R = t_y/t_x$ is linearised as $u = (y - \hat R x)/\hat t_x$ (Woodruff 1971), and a mean is the ratio with $x = 1$.

## 5. Survival and Censoring

*Module `survival`.*

**Theorem 5.1 (Kaplan and Meier 1958).** Under independent censoring, the hazard of the event at an event time $t_j$ among those at risk is estimated by $d_j/n_j$, and $S(t) = \prod_{t_j \le t} (1 - d_j/n_j)$ is the nonparametric maximum likelihood estimator of the survival function. Its asymptotic variance is estimated by Greenwood's formula $S^2 \sum d_j/(n_j(n_j - d_j))$.

*Remark (censoring as coarsening).* With $G(t) = P(C > t)$, the Kaplan–Meier estimate of $1 - S(t)$ equals the inverse-probability-of-censoring weighted estimate $\frac1n \sum_i \delta_i [X_i \le t]/\hat G(X_i-)$, with $\hat G$ the Kaplan–Meier estimate of the censoring law (Satten and Datta 2001). It is the Horvitz–Thompson estimator of §2, with being uncensored at one's event time as the inclusion. This identity is planned as a check.

**Theorem 5.2 (Cox 1972; Efron 1977).** Under $\lambda(t \mid x) = \lambda_0(t) e^{x\beta}$, the partial likelihood $\prod_j e^{x_{(j)}\beta} / \sum_{i \in R(t_j)} e^{x_i\beta}$ does not depend on $\lambda_0$. With ties, Breslow's approximation uses the full risk set for every tied event. Efron's removes the tied events' share $l/d$ step by step. The score test at $\beta = 0$ of one binary covariate, with no ties, is the log-rank statistic.

*Checked* against statsmodels and lifelines on Rossi's study to $10^{-13}$. The score/log-rank identity is checked exactly.

**Theorem 5.3 (Aalen and Johansen 1978).** With competing causes, the cumulative incidence $F_c(t) = \sum_{t_j \le t} S(t_j -) d_{cj}/n_j$ satisfies $\sum_c F_c(t) = 1 - S(t)$.

*Proof.* $1 - S(t_j) = (1 - S(t_j-)) + S(t_j-) d_j/n_j$ with $d_j = \sum_c d_{cj}$. Induct on $j$. ∎ *Checked* to $10^{-16}$.

## 6. Missing Data and Multiple Imputation

*Module `missing`.*

**Definition 6.1 (MAR).** Nonresponse is *missing at random* when $P(R = 1 \mid Y, X) = P(R = 1 \mid X)$, CAR for the nonresponse row of §1. Under MAR and separate parameters the response mechanism is *ignorable* for likelihood inference (Rubin 1976).

**Theorem 6.2 (Rubin's rules; Rubin 1987).** Let $m$ imputations be drawn properly from the posterior predictive law of the missing values. Let $\hat Q_l$, $U_l$ be the complete-data estimates and variances. Then $\bar Q = \frac1m \sum \hat Q_l$ is the point estimate, with total variance $T = \bar U + (1 + 1/m) B$, $B$ the between-imputation variance, and $(\bar Q - Q)/\sqrt{T}$ referred to a $t$ law on $\nu = (m - 1)/\gamma^2$ degrees of freedom, $\gamma = (1 + 1/m) B / T$. Barnard and Rubin (1999) correct $\nu$ for complete data with few degrees of freedom.

The approximate Bayesian bootstrap (Rubin and Schenker 1986) is a proper imputation within classes: resampling the donors first carries the uncertainty of each class's law. *Checked* by simulation: under MAR given a class, the complete-case mean is biased by $-1.16$, while twenty imputations pooled by Rubin's rules cover the mean at the nominal rate.

## 7. Identification by Assumptions on the Outcome

*Module `quasi`.* When the mechanism of assignment is not known, an assumption on the potential outcomes can take its place. Each result below has the shape of Definition 1.2, with an assumption that the data cannot test entirely.

**Theorem 7.1 (group-time effects; Callaway and Sant'Anna 2021, Thm 1).** Let units be first treated in period $G \in \{2, \dots, T\} \cup \{\infty\}$, treatment being absorbing, and let $C$ be a comparison group (never treated, or not yet treated by $t$). Suppose there is no anticipation ($Y_t = Y_t(0)$ for $t < G$), and parallel trends: $E[Y_t(0) - Y_{g-1}(0) \mid G = g] = E[Y_t(0) - Y_{g-1}(0) \mid C]$. Then for $t \ge g$
$$\mathrm{ATT}(g, t) = E[Y_t(g) - Y_t(0) \mid G = g] = E[Y_t - Y_{g-1} \mid G = g] - E[Y_t - Y_{g-1} \mid C].$$

*Proof.* $E[Y_t - Y_{g-1} \mid G = g] = \mathrm{ATT}(g, t) + E[Y_t(0) - Y_{g-1}(0) \mid G = g]$ by no anticipation. Substitute parallel trends, and use that the comparison units are untreated at $t$ and $g - 1$. ∎

*Remark.* The static two-way regression is a weighted sum of such comparisons, some with negative weights when effects vary with exposure (Goodman-Bacon 2021). *Checked:* with effects growing with exposure in two groups, the regression gives 1.94 where the average effect is 2.43, while every $\mathrm{ATT}(g, t)$ is recovered exactly. With one constant effect the two agree.

**Theorem 7.2 (the local average treatment effect; Imbens and Angrist 1994).** Let $Z \in \{0, 1\}$ be an instrument and $D(z)$, $Y(d)$ the potential treatments and outcomes, with $Z$ independent of $(D(0), D(1), Y(0), Y(1))$. Suppose $Z$ affects $Y$ only through $D$, monotonicity holds ($D(1) \ge D(0)$), and $P(D(1) > D(0)) > 0$. Then
$$\frac{E[Y \mid Z = 1] - E[Y \mid Z = 0]}{E[D \mid Z = 1] - E[D \mid Z = 0]} = E[Y(1) - Y(0) \mid D(1) > D(0)].$$

*Proof.* The numerator is $E[(D(1) - D(0))(Y(1) - Y(0))]$, and with monotonicity $D(1) - D(0)$ is the indicator of compliers. The denominator is their share. ∎

*Checked:* the 2SLS estimate is the ratio exactly, and recovers the compliers' effect beside always-takers whose effect differs.

**Theorem 7.3 (regression discontinuity; Hahn, Todd and van der Klaauw 2001).**
- **Sharp.** If $D = [x \ge c]$ and $E[Y(0) \mid x]$, $E[Y(1) \mid x]$ are continuous at $c$, then $\lim_{x \downarrow c} E[Y \mid x] - \lim_{x \uparrow c} E[Y \mid x] = E[Y(1) - Y(0) \mid x = c]$.
- **Fuzzy.** If $P(D = 1 \mid x)$ jumps at $c$, the outcomes' means are continuous, and monotonicity holds near $c$, then the ratio of the jumps in $E[Y \mid x]$ and $E[D \mid x]$ is the compliers' average effect at $c$.

Local linear regression estimates each one-sided limit with bias $O(h^2)$ at the boundary, without the $O(h)$ of a local mean (Fan and Gijbels 1996). *Checked:* the bias falls by 3.9998 when $h$ halves, lines give the jump exactly, and the fuzzy estimate is the ratio of the jumps.

**Synthetic control.** When outcomes follow a linear factor model $Y_{it}(0) = \delta_t + \lambda_t \mu_i + \varepsilon_{it}$, the synthetic control's bias is bounded by terms that vanish as the pre-treatment fit becomes exact and $T_0$ grows (Abadie, Diamond and Hainmueller 2010, §2.3). Inference is by permutation in space. If the treated unit was as likely as any donor to be the one treated, the rank of its post-to-pre ratio has the uniform law on $\{1, \dots, J + 1\}$, so the p-value is exact for that null (ibid., §2.4).

## 8. Causal Graphs

*Module `graphs`.* A directed acyclic graph $\mathcal{G}$ on variables encodes a causal model, $X_v = f_v(\mathrm{pa}(v), \varepsilon_v)$ with independent $\varepsilon$.

**Definition 8.1 (d-separation; Pearl 1988).** A path is *blocked* by $Z$ when it contains a non-collider in $Z$, or a collider that is not in $Z$ and has no descendant in $Z$. $X$ and $Y$ are *d-separated* by $Z$ when every path between them is blocked.

**Theorem 8.2 (soundness and completeness).**
- **Soundness.** If $X$ and $Y$ are d-separated by $Z$, then $X \perp Y \mid Z$ in every distribution that factorises over $\mathcal{G}$ (Verma and Pearl 1988).
- **Completeness.** If they are not, the independence fails in some such distribution (Geiger, Verma and Pearl 1990). In linear Gaussian models it fails for all coefficients outside a set of measure zero (Spirtes, Glymour and Scheines 2000).

*Checked:* on twelve random graphs with random coefficients, d-separation and the vanishing of the partial correlation computed from $\Sigma = (I - B)^{-1}(I - B)^{-\top}$ agree in all 4,032 cases.

**Theorem 8.3 (adjustment; Pearl 2009, Thms 3.3.2 and 3.3.4).**
- **Back door.** If $Z$ contains no descendant of $X$ and d-separates $X$ from $Y$ in $\mathcal{G}$ with the arrows out of $X$ removed, then $P(y \mid do(x)) = \sum_z P(y \mid x, z) P(z)$.
- **Front door.** If $Z$ intercepts every directed path from $X$ to $Y$, no back-door path from $X$ to $Z$ is open, and $X$ blocks every back-door path from $Z$ to $Y$, then $P(y \mid do(x)) = \sum_z P(z \mid x) \sum_{x'} P(y \mid x', z) P(x')$.

*Checked* on the elementary graphs, on M-bias (adjusting for the collider $M$ opens a path) and on the front door with a hidden confounder.

## 9. What Comes Next

These results are written in the later stages of the infrastructure's ADR 0082:
- **Cohorts** (begun, the module `cohorts`): the Lexis diagram's person-time and events, exact rates, standardisation, for demography (the infrastructure's ADR 0085). Next in it is non-identification shown exactly: in the age–period–cohort model, age $=$ period $-$ cohort, so the design matrix loses one dimension of its linear trends. Only the second differences (curvatures) are identified (Holford 1983).
- **Partial identification.** Without assumptions on the mechanism, effects are bounded, not identified (Manski 1990). The bounds and Rosenbaum's sensitivity parameter quantify how much an untestable assumption carries.
- **Synthesis and sequential testing:** meta-analysis, and boundaries that keep their level when data are looked at as they come.
