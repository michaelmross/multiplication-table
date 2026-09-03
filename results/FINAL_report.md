# The per-gap multiplication-table density from N = 10³ to N = 2^(10^12): final report

## Campaign structure

Three instruments measured one quantity — ρ(N, c), the density of integers near (cN)² possessing a divisor within factor 1/c of their square root, equivalently the distinct-entry density of the N × N multiplication table in the square-gap at k = cN. The exact census (divisor sieve, N ≤ 10⁸, brute-force verified) covers t = log log x up to 3.6. The Tier-1 exact-integer Monte Carlo (Kalai random factored integers, gmpy2, run on the author's hardware across multiple interrupted-and-resumed sessions including one filesystem fault) extends exact arithmetic to N = 2^1024, t = 7.26. The Tier-2 continuum model (GEM(1) stick-breaking factorizations, the uniform-order-statistics anatomy underlying Ford's theorem) runs to N = 2^(10^12), t = 27.96, anchored to the exact tiers at four seam heights. Predictions P1–P5 were preregistered in the script headers with evaluation criteria frozen in `analyze_mc.py` before any Tier-1 or final Tier-2 data existed.

## The exact-integer ladder (Tier 1, c = 1/2)

| height | ln x | ρ | rel. SE | samples | Kalai attempts / (2.82 ln x) |
|---|---|---|---|---|---|
| N = 10⁸ | 35.5 | 0.24973 | 1.3% | 17,875 | 1.001 |
| 2^64 | 87.3 | 0.19979 | 3.3% | 3,779 | 0.974 |
| 2^128 | 176.1 | 0.17119 | 4.0% | 2,985 | 0.986 |
| 2^256 | 353.5 | 0.14448 | 4.4% | 2,997 | 0.993 |
| 2^512 | 708.4 | 0.12039 | 4.6% | 3,389 | 1.025 |
| 2^1024 | 1418.2 | 0.11407 | 12.2% | 526 | ~0.94 |

The attempts-per-sample column is a continuous health certificate: Kalai's acceptance probability is Mertens' product, so attempts per kept sample must equal e^γ/(1 − e^(−1)) · ln x = 2.818 ln x, and the data sits on that constant within ±3% across a 40-fold range of ln x. The same e^(−γ) that appeared as the 2e^(−γ) Mertens overshoot in the R_N census reappears here as the price of exact sampling.

**P1 (harness validity): PASS.** Tier-1 at N = 10⁸ gives 0.24973 against the exact census value 0.25276 — deviation 0.0030 against tolerance 0.0065.

## The seams (model vs exact arithmetic)

At c = 1/2, seven independent comparisons of the continuum model against exact truth (census values at 10⁴, 10⁶, 10⁸; Tier-1 at 2^64, 2^256, 2^512, 2^1024) give offsets −8.1%, −5.2%, −0.2%, −7.2%, −8.0%, −0.1%, −9.4%, with seam standard errors of 6–15 percentage points on the Tier-1 entries. The offsets scatter about a constant ≈ −5% with no detectable height trend across four decades of ln x. This is the load-bearing empirical fact of the construction: a height-stable multiplicative offset cancels identically in d(ln ρ)/dt, so the Tier-2 slope results below carry exact-integer authority.

**P2 (model offset ≤ 20% everywhere): FAIL as stated.** The bound holds at c = 0.5 but the narrow windows breach it: the model runs −19% to −23% below exact truth at c = 0.75 and c = 0.9, height-stable. The diagnosis stands as recorded at preregistration time: the continuum model carries no prime multiplicities and no small-prime discreteness, and the missing divisors cost most where the window is narrowest. The failure is a measured property of the model, localized, and — because it is height-stable — harmless to the slope program.

## Shape emergence (P3, P4)

**Both FAIL as preregistered, and the failures measure something.** The preregistered threshold said Ford's local shape −δ − (3/2)/t would hold within 0.05 from t = 8. In fact the deviation — which begins at +0.36 in the exact-census regime — decays through +0.10 around t ≈ 8–11 and enters the 0.05 band durably only at t ≈ 13–16, roughly N of a few million bits. From there to the top of the range the shape holds at every c, with two excursions of 0.054 and 0.051 that are individually within 2σ. Likewise c-universality (P4): the preregistered pairwise band of 0.03 is met on t ∈ [16, 25] (spreads 0.014–0.025) but not below, and the single breach above (spread 0.038 at t = 26.8, with per-point SEs of 0.02) is noise-compatible. The honest one-line summary: the Ford asymptotic shape is numerically real and universal, and it turns on three e-folds of log log later than guessed.

## The exponent (P5)

**PASS by the letter, and the letter now conceals the finding.** The pooled fit over t ∈ [10, 28] gives d̂ = 0.0596 ± 0.0042 against δ = 0.086071, inside the preregistered ±0.03 band — but the author's 2500-hit rerun reduced the noise so far below the anticipated level that the band and the error bar decoupled, and the point estimate is six standard errors below δ. Windowing resolves what happened:

| t-window | d̂ (pooled, all c) |
|---|---|
| [4, 10) | −0.048 ± 0.011 |
| [10, 15) | +0.031 ± 0.014 |
| [15, 20) | +0.042 ± 0.007 |
| [20, 28.5) | +0.074 ± 0.006 |

The effective exponent climbs monotonically toward δ and has not arrived: even in the top window — heights 2^(10^8) to 2^(10^12), trillion-bit integers — it sits 2.2σ short. Since the (3/2)/t term is confirmed by the shape tests, the residual is a genuine subleading correction to Ford's asymptotic in the model, still visible at 10^12 bits and decaying slowly enough to be consistent with an O(1/t)-type term of moderate coefficient. Ford's theorem is an order-of-magnitude statement and asserts nothing about this structure, so the observation is a numerical characterization of unexplored territory: **the multiplication-table exponent δ is not attained, even by its own generating mechanism, at any height ever likely to be touched by exact computation — its asymptopia lies beyond N = 2^(10^12).**

## Reading for the program

The campaign closes the question that opened it — whether individual square-gaps track Ford's average, and where the asymptotics live. The answer has three layers. Gap-to-gap, the density is rigid with square-root cancellation (CV ~ N^(−1/2), the strongest confirmed prediction of the exact census), so the unproven localization of Ford's theorem to √x-scale windows is empirically a concentration statement holding with room to spare. Level-wise, the density's decay follows the Ford shape only from millions of bits upward, and its exponent remains measurably shy of δ at a trillion bits. Methodologically, the whole structure rests on two constants of Mertens type showing up exactly where theory puts them — 2e^(−γ) in the census, e^γ ln x in the sampler — which is what licenses trusting the instruments about everything else.

For the Jₙ program the moral is the one suspected from the start, now with numbers attached: at √x scale the arithmetic fluctuations are tame and the obstructions are entirely in the methods, while the constants that govern the levels converge so slowly that no computational height can substitute for a proof. The census-to-model pipeline built here — exact sieve, exact-integer sampler, anchored continuum model, preregistered throughout, with three refutations and two passes honestly logged — is reusable for any divisor-anatomy quantity in the program, and the subleading-correction observation at P5 is, as far as I can tell, new.

## Provenance

Exact census: `gap_census.py`, `gap_summary.csv` (brute-force verified). Tier 1: `t1_exact_mc.py`, `mc_tier1.csv`, 93 rows across resumed sessions, 32,551 kept samples, seeds advancing per round. Tier 2: `t2_model_mc.py`, `mc_tier2.csv`, 48 cells, ~600–2500 hits per cell. Analysis: `analyze_mc.py` (criteria frozen; pooling of resumed rounds by sample count). Charts: `mc_chart_slopes.png` (local decay rate against Ford shape across 24 e-folds of t, all three instruments), `mc_chart_density.png` (the density itself). Design amendments after profiling but before measurement are logged in the script headers, including one corrected model bug (window centering) and one corrected overflow bug (integer keep-cut past the float64 ceiling, first triggerable at exactly the 2^1024 height).


## Erratum and refinement (post-report, logged before Campaign A)

Preparing the deep-t extension exposed that the three overlap anchors (N = 10^4, 10^6, 10^8) ran only their 3,000-sample pilots in the August campaign, because at overlap densities the 600-hit goal sits below the pilot size. Their offsets therefore carried ±2.5-2.7 point errors that the report's narrative treated as sharper than they were. High-n top-ups (80,000 samples each, c = 0.5, engine verified equivalent to the original on 40,000 identical inputs with zero disagreements) refine the c = 0.5 model-vs-exact offsets to a height-stable deficit of about -5 to -7% with no residual suggestion of shrinkage toward zero at N = 10^8; the earlier "-0.2% at 10^8" was pilot noise. This strengthens rather than weakens the report's load-bearing claim: the c = 0.5 offset is constant in height, so Tier-2 slopes inherit exact-integer validity. No verdict changes. The top-up rows are appended to mc_tier2.csv with n = 80,000 and pool correctly in analyze_mc.py.


## Campaign B results (ensemble exact census, author's hardware)

Both preregistered criteria pass. B1 (square-centered typicality, now
powered): the TV distance of J_n's multiplicity histogram from the 32-window
shifted-ensemble mean sits at z = -0.02, +0.90, -0.26 for N = 10^6, 10^7,
10^8. J_n is statistically indistinguishable from a generic window of its
magnitude at every anchor, retiring the single-window caveat of the original
census with a genuine null test. B2 (concentration exponent): CV(D) over 50
consecutive square-gaps fits CV ~ N^beta with beta = -0.440 +- 0.031
(computed from the reported summary values at two significant figures,
pending the full-precision CSV; the rounding perturbs beta by at most
~0.015). The estimate lies inside the preregistered band [-0.6, -0.4], is
+1.9 sigma from exactly -1/2, and the per-decade slopes
(-0.469, -0.410) are
individually in band. Combined verdict: the unproven localization of Ford's
theorem to sqrt(x)-scale windows is, empirically, a concentration statement
with a fitted square-root exponent and a powered typicality null behind it.


## Campaign A results (deep-t extension to N = 2^(10^15), author's hardware)

Frozen verdicts: A1 met (pooled deficit on t in [28, 35.5] = +0.0137 +/-
0.0032, SE target satisfied). A2: z0 = +4.29 against H_0 -- **the deficit is
real; the continuum model has not reached its asymptotic exponent at
quadrillion-bit heights** -- and z1 = +1.52 against the frozen 1/t
hypothesis, so H_1/t survives. Joint fit over all t >= 20: b = 0.359 +/-
0.084. The effective exponent at t ~ 31 is delta - 0.0137 = 0.0724 +/-
0.0032, i.e. 4.3 sigma below delta = 0.086071.

Post-hoc caution, logged as such: a constant-deficit alternative (not among
the frozen hypotheses) fits the two t >= 20 windows as well as 1/t does --
the pooled deficit moved +0.0122 -> +0.0137 between window midpoints 24 and
31, flat within errors. The three-window picture (deficit 0.044 at t ~ 17.5,
then ~0.013 flat) rules out a global constant and a global 1/t alike, and is
consistent with a faster-decaying term dominating below t ~ 20 plus either a
slowly-decaying tail or a constant offset above it. A constant offset would
mean the model's asymptotic exponent differs from Ford's delta -- surprising
if the GEM mechanism is truly the delta-generating one, but that
universality is a belief, not a measurement, and the data currently cannot
reject it. Discriminating b/t from constant at 3 sigma requires either
roughly quadrupled hits on the t in [20, 28.5] cells (2^10^9 .. 2^10^12; the
cheaper lever, since that window's SE dominates the comparison) or one more
lever arm at t ~ 40, which exceeds longdouble headroom and would need the
integer+fraction exact-window representation. Either is a preregisterable
Campaign A-prime.


## Campaign A' results and instrument erratum (N to 2^(10^18))

The frozen D2 verdict returned "excludes BOTH hypotheses" -- and forensics
attribute this to an instrument failure, not arithmetic of divisors. The
anomaly concentrates in the last interval and is far larger at the narrow
window (c = 0.75 deficit +0.0587 at t in [39.5, 41.8]), the signature of
precision blur in the window test. The generator was exonerated by a
theory-anchored diagnostic (the stick intensity dl/l fixes the expected
count in [ln 2, e ln 2) at exactly 1 at every height; additive and
multiplicative remainder tracking agree to within SE at all L). A blur
dose-response test at a clean height then showed that O(0.2) noise in the
window comparison -- the longdouble rounding scale at L ~ 10^18 -- moves the
measured density by tens of percent. The sparsity argument that certified
longdouble to 2^(10^18) was therefore wrong: it ignored the clustering of
subset sums around near-hits. The 2^(10^17) and 2^(10^18) cells are voided
pending an exact-arithmetic rerun; the author's compute on them is the cost
of the instrument lesson.

On arithmetically clean heights (t <= 38): the deep point at t_bar = 36.0
gives deficit +0.0070 +- 0.0027 -- z = -2.10 against the
frozen constant floor, z = -1.10 against the frozen 1/t --
and together with the 4.8-sigma decline between the mid windows, the
constant-floor hypothesis is effectively dead while 1/t survives. The clean
two-parameter fit gives C_inf = -0.0092 +- 0.0035
(consistent with zero) and b = +0.623 +- 0.085;
the pure 1/t fit gives b = +0.403. Provisional physical conclusion,
pending the exact-engine rerun of the top arm: the continuum model's
asymptotic exponent is Ford's delta, approached through a subleading
correction consistent with b/t, b ~ 0.3.


## Arbiter verdict: blur artifact confirmed with quantitative closure

The 60,000-pair exact-vs-longdouble arbiter returned seven classification
flips across the deep heights -- all seven longdouble-only (p ~ 0.016 for
directionality by sign test): longdouble blur widens the window and INFLATES
the measured density, by +1-2% relative at 2^10^18 (c = 0.5) and +15.5% at
2^10^18 (c = 0.75), with 2^10^17 nearly clean. Propagating this inflation
through the slope arithmetic predicts a spurious deficit excess on the
[39.5, 41.8] interval of about +0.009 (c = 0.5) and +0.052 (c = 0.75);
the observed anomalous excesses were +0.011 and +0.050. The artifact
explains the voided cells to within the flip-count Poisson errors, in both
magnitude and c-dependence. The "real physics" reading of the deep-arm
anomaly is closed; the earlier jitter dose-response had in fact shown the
correct (inflationary) sign. The exact engine's gate passed at 100.0000%
agreement on the clean height; the rerun proceeds on a validated instrument.


## Campaign A3 final verdicts (exact engine, N to 2^(10^18))

The exact-arithmetic deep arm landed on the frozen 1/t prediction to within
measurement resolution: deficit +0.0090 +- 0.0022 at t_bar = 39.5 against
the frozen prediction +0.0091 -- z = -0.06. Frozen D2 as written: H_1/t
consistent; H_const disfavored at the deep arm alone (z = -1.68). The joint
comparison over all 18 clean per-interval points with t >= 20 settles
it: the best constant floor pays chi2 = 88.4 against 31.6 for
b/t (each with 17 dof), delta-chi2 = 56.8, rejecting the constant
at ~7.5 sigma equivalent, while b/t fits with
chi2/dof = 1.86. The coefficient is stable across the mature
regime: window values b = 0.344 +- 0.038 (t = 31.3) and 0.355 +- 0.087
(t = 39.5); global fit b = +0.402 +- 0.014; restricted to t >= 28,
b = +0.333 +- 0.032 with chi2/dof = 1.12. The
two-parameter fit's mildly negative C_inf (-0.0074 +- 0.0030) is the
sub-t=25 transient bending a two-parameter shape; on mature windows C_inf
is consistent with zero.

**Physical conclusion of the program.** The continuum (GEM) model of the
multiplication-table density approaches Ford's exponent delta = 0.086071 as
its true asymptotic exponent -- no floor -- through a subleading correction
consistent with b/t, b ~ 0.35, preceded by a faster transient dying around
t ~ 25 (N ~ 2^(10^10)). The local shape -delta - (3/2)/t is confirmed from
t ~ 13-16. Delta itself is never attained at any computed height: at
t = 39.5 (N = 2^(10^18), quintillion-bit integers) the effective exponent
still reads 0.077. All of it measured on a validated instrument, with the
voided longdouble cells explained to Poisson precision by the paired
arbiter, and every frozen criterion evaluated as written.

## Companion computation: the ballot skeleton (pre-submission check)

Prompted by the question of whether the slow approach was mechanism-forced
and cheaply computable: largely yes, and the agreement is the finding. The
exactly computable ballot skeleton -- a unit-rate Poisson process on [0, t]
required to stay above the line of slope 1/ln2, i.e. Ford's uniform-order-
statistics mechanism with all divisor structure stripped away -- is a
50-line positive-arithmetic DP validated against Monte Carlo. It reproduces
delta as a large-deviation rate, the universal t^(-3/2), and a subleading
deficit of +0.0091..0.0109 at t ~ 31 and +0.0063..0.0072 at t ~ 40 (theta
in [0.5, 2] as boundary systematic), against the campaign's measured
+0.0110 +- 0.0012 and +0.0090 +- 0.0022; its effective 1/t coefficient on
the measured range is 0.26-0.31 against b = 0.333 +- 0.032. Conclusion: the
1/t correction class and (within ~1-2 sigma) its coefficient are properties
of the ballot mechanism itself, and the lossy steps in the reduction from
divisors to the ballot event do not carry the subleading term. Where model
and skeleton genuinely differ is the transient (model ~0.040 vs skeleton
~0.025-0.029 at t ~ 16): the divisor structure the skeleton drops lives in
the fast-decaying part, now localized. The skeleton's coefficient slides to
0.21-0.24 on t in [45, 62], yielding a falsifiable prediction -- b should
decline to ~0.22 by t ~ 50 -- for any future campaign beyond the current
exact-engine ceiling. The note's framing changes accordingly: the campaign
did not discover that the approach is slow (mechanism-predicted, computable
in seconds); it established that the full model transfers the skeleton's
subleading law, anchored it to exact integers at the seams, measured the
coefficient, and localized the transfer's failure to the transient.
