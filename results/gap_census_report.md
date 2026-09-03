# Distinct table entries per square-gap: preregistered results

**Verdict up front: P4 confirmed strongly; P1, P2, P3 refuted as stated.** The refutations are informative and diagnosable, and the diagnosis is corroborated by the published numerics on the multiplication-table problem. Details below; nothing is reframed post hoc — the predictions are quoted as preregistered in `gap_census.py`, which also contains the brute-force correctness cross-check (exact table construction at N = 97, 100, 211, 300, all c; passed before the main run).

**Object.** For k = ⌊cN⌋ and the square-gap G_k = (k², (k+1)²], D(N, c) counts distinct entries of the N × N multiplication table in G_k, via the membership criterion m ∈ table ⟺ m has a divisor in [m/N, N]; ρ(N, c) = D/(2k+1). Ensembles of 10 consecutive gaps (3 at N = 10⁸) give per-gap dispersion. Anchors N = 10³ … 10⁸; c ∈ {0.5, 0.75, 0.9}; exact integer arithmetic in the sieve.

## Results

| N | ρ (c=0.5) | ρ (c=0.75) | ρ (c=0.9) | CV (c=0.5) | ρ/ρ₀ (c=0.5) | Ford factor |
|---|---|---|---|---|---|---|
| 10³ | 0.34060 | 0.21544 | 0.09448 | 1.94% | 1.0000 | 1.0000 |
| 10⁴ | 0.31163 | 0.20239 | 0.09285 | 0.38% | 0.9149 | 0.7922 |
| 10⁵ | 0.29056 | 0.19110 | 0.08971 | 0.20% | 0.8531 | 0.6732 |
| 10⁶ | 0.27495 | 0.18258 | 0.08700 | 0.07% | 0.8073 | 0.5949 |
| 10⁷ | 0.26270 | 0.17577 | 0.08472 | 0.02% | 0.7713 | 0.5389 |
| 10⁸ | 0.25276 | 0.17017 | 0.08280 | 0.01% | 0.7421 | 0.4965 |

**P1 (Ford-shape decay): refuted.** The preregistered claim was agreement of ρ(N)/ρ(10³) with F(N) = (log N/log 10³)^(−δ)(log log N/log log 10³)^(−3/2) within 20% at every anchor. Observed at N = 10⁸: 0.742 (c = 0.5), 0.790 (c = 0.75), 0.876 (c = 0.9) against F = 0.497 — deviations of 49%, 59%, 76%. The empirical decay is real (the constant-density null, predicting 1.00, also fails) but far slower than the asymptotic shape.

**P2 (universality in c): refuted.** The normalized decay curves do not coincide: the far-anchor values 0.742/0.790/0.876 separate monotonically in c, i.e., the c-dependence has not factored out of the N-dependence at these heights. The level ordering (ρ decreasing in c, wider window ⟹ denser) holds as predicted, but that was the weak half of P2.

**P3 (slope diagnostic): refuted.** Fitted d with the (log log)^(3/2) factor imposed came out at −0.32, −0.38, −0.49 — negative, meaning that after dividing out (log log N)^(−3/2) the residual *grows* with log N. The imposed 3/2 factor overshoots the entire observed decay at these heights.

## Diagnosis

The three refutations share one cause, visible in the local-slope chart: in the variable t = log log N, Ford's shape predicts a local decay rate d log ρ/dt = −δ − (3/2)/t ≈ −0.61 over our range, while the empirical rates on the last interval are −0.289, −0.243, −0.172 (c = 0.5, 0.75, 0.9) — between a quarter and a half of the asymptotic rate, and still c-dependent. What is refuted is not Ford's theorem (an order-of-magnitude statement with unspecified constants, proven) but the auxiliary stationarity assumption embedded in P1–P3: that the slowly varying factor multiplying the asymptotic shape is already constant for N ∈ [10³, 10⁸]. It is not; the effective constant is still climbing throughout the accessible range, absorbing most of the (log log)^(3/2) decay. At our anchors t = log log N only moves from 1.93 to 2.91 — asymptotically, "t → ∞" has barely begun.

This is corroborated independently: the Brent–Pomerance–Purdum–Webster computations of M(n) ran exact values to 2³⁰ and Monte Carlo estimates to heights as large as 2^(10⁸) precisely because the asymptotic behavior of M(n) is not resolvable at exact-computation heights, and even their comparison against Ford's order-of-magnitude bound is a statement about slowly drifting normalized ratios, not a measured limit (arXiv:1908.04251). Our per-gap census reproduces the same glacial pre-asymptotics in a localized setting. A follow-up worth considering: a Monte Carlo per-gap analogue at Brent–Pomerance heights (sampling m ∈ G_k and testing for a divisor in [m/N, N] via Bach/Kalai-style random factored integers), which would test whether the local slope bends toward −δ − (3/2)/t as t grows — the falsifiable version of "the pre-asymptotic drift is the whole story."

**P4 (per-gap concentration): confirmed, strongly.** The coefficient of variation of D over consecutive gaps falls from 1.94% at N = 10³ to ≤ 0.01% at N = 10⁸, tracking N^(−1/2). Individual square-gaps hug their local mean to a few parts in 10⁴ at the top anchor. Empirically, the "almost all gaps" statement that no one has proved is not merely true but true with square-root cancellation — consecutive gaps at scale k ≍ cN are, to this resolution, exchangeable samples of one local density. The unproven localization of Ford's theorem is, on this evidence, a concentration statement waiting for a second-moment method strong enough to reach length-√x windows; the obstruction to proving it is the usual one, but the census now says the target is not merely plausible but tight.

## Combined reading with the R_N census

The two censuses now bracket the phenomenon from both sides. The multiplicity census showed the *strip* statistics of J_n are exactly generic (TV at noise floor, PNT rate at R = 1, the 2e^(−γ) Mertens overshoot in its textbook place). The gap census shows the *distinct-entry* statistics of square-gaps are locally rigid (CV ~ N^(−1/2)) while their absolute level drifts through a pre-asymptotic regime that five decades cannot exit. Jointly: localization at √x scale looks empirically like a concentration problem with strong square-root cancellation and a slowly varying local density — the fluctuation side is tame; it is the methods, not the arithmetic, that cannot yet certify it.

## Files

`gap_census.py` — preregistered script with brute-force cross-check. `gap_summary.csv` — all anchors, all c, ensembles. `gap_chart_density.png` — ρ vs N with Ford-shape overlays (divergence visible). `gap_chart_slope.png` — local decay rate vs Ford's −δ − (3/2)/t (the diagnosis chart). `gap_chart_cv.png` — P4 concentration with N^(−1/2) guide.

Runtime: 27 s total for all exact counts including ensembles; peak memory ~180 MB.
