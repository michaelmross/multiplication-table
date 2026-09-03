# Strip-multiplicity census for J_n: preregistered results

**Object.** For J_n = [4n² − n, 4n² + n] and N = 2n, define the strip multiplicity R_N(m) = #{a ≤ N : a | m}, the number of appearances of m in the N-row strip of the multiplication table (counting row 1). Since every composite in J_n has a prime factor ≤ N and no element of J_n factors as a product of two integers both exceeding N, we have R_N(m) = 1 if and only if m is prime. Legendre near squares is therefore a statement about the lower tail of the multiplicity distribution of the multiplication table restricted to a width-2n window at the corner.

**Design.** Four predictions were stated in the script header before execution (P1–P4 in `rn_census.py`). Anchors n = 10³, 10⁴, 10⁵, 10⁶, 10⁷; the largest window has width 2·10⁷ + 1 centered at 4·10¹⁴. The comparison window S_n = [4n² + 10n, 4n² + 12n] is a generic shifted window of identical width at the same magnitude. R was computed by direct divisor sieve (exact integer arithmetic; the only floating point is in the reported ratios and the Mertens product, accumulated via log1p).

## Preregistered results

| n | primes (R=1) | #{R=1}/u′ (P2) | Mertens/emp (P3) | mean R − H_N (P1) | TV(J_n, shifted) (P4) |
|---|---|---|---|---|---|
| 10³ | 131 | 0.99522 | 1.12532 | +0.01903 | 0.02599 |
| 10⁴ | 1,019 | 1.00911 | 1.11217 | +0.00495 | 0.01600 |
| 10⁵ | 8,162 | 0.99625 | 1.12693 | −0.00112 | 0.00482 |
| 10⁶ | 68,916 | 0.99988 | 1.12300 | +0.00057 | 0.00143 |
| 10⁷ | 595,645 | 1.00135 | 1.12139 | −0.00003 | 0.00055 |

Target for P3: 2e^(−γ) = 1.12292.

**P1 (mean).** Confirmed, and more sharply than predicted: the deviation mean R − H_N is not merely O(1) but decays roughly like 1/√n, indicating that the per-row boundary errors O(1) in ⌊·⌋-counting are equidistributed in sign across rows and cancel on average. The mean of the multiplicity distribution is pinned to H_N = log N + γ + o(1) with no visible bias from square-centering.

**P2 (PNT rate).** Confirmed. The prime count in J_n sits at the PNT rate u′ = (2n+1)/log(4n²) within ±1% at every anchor and within 0.14% at the two largest.

**P3 (Mertens overshoot).** Confirmed. The independent-divisibility model — treating divisibility by each prime p ≤ N as an independent event of probability 1/p — predicts (2n+1)·∏_{p≤N}(1 − 1/p) integers with R_N = 1, and this overshoots the truth by a factor converging to 2e^(−γ): the ratio reads 1.12300 at n = 10⁶ against the target 1.12292, agreeing to four decimal places, with the residual oscillation at 10⁷ consistent with the ±0.1% sampling noise visible in P2. This is the numerical face of the fact that the sieve's independence heuristic fails by a universal constant precisely at the R = 1 level — the toy shadow, in this window, of why moment and independence methods cannot certify the lower tail and must instead be supplied with equidistribution input.

**P4 (window typicality).** H₀ retained. The total-variation distance between the multiplicity distributions of the square-centered and shifted windows decays like n^(−1/2), i.e., exactly at the sampling-noise floor for windows of this width. At the resolution of the full R_N-distribution, square-centering is invisible: whatever is special about J_n (the divisor pairing near √m, the self-conjugate parabola geometry) does not register in first-order divisor statistics. This is consistent with, and a distributional restatement of, the coefficient-frequency finding in the Jₙ program: the obstruction is not visible in unsigned counts.

*Caveat recorded in advance:* R_N = 1 on the shifted window means "prime or p² with p ∈ (2n, 2n+3]," a discrepancy of size O(1) per window, and the shifted comparison uses a single window per anchor, so P4 is a diagnostic against gross anomalies, not a powered test. A follow-up with an ensemble of shifted windows would put error bars on the TV floor.

## Exploratory (not preregistered)

**R = 2 stratum.** The count of m ∈ J_n with R_N(m) = 2 — up to O(1) prime-power cases, the semiprimes pq with p ≤ N < q, i.e., the P₂ stratum that Chen-type methods reach and the stratum whose tail T(θ) the one-hypothesis paper bounds — tracks the heuristic u′·(log log N + M + log 2) with a slowly decaying deficit of expected 1/log-type size:

| n | R=2 count | (R=2)/u′ | log log 2n + M + log 2 |
|---|---|---|---|
| 10³ | 375 | 2.8489 | 2.9829 |
| 10⁴ | 3,180 | 3.1492 | 3.2475 |
| 10⁵ | 27,755 | 3.3878 | 3.4566 |
| 10⁶ | 246,572 | 3.5774 | 3.6294 |
| 10⁷ | 2,221,054 | 3.7339 | 3.7767 |

**Maximum multiplicity.** The champions max R_N grow from 96 (n = 10³) to 5,376 (n = 10⁷); doubled, these are divisor counts 192 → 10,752 of the most composite integers in each window, the upper tail complementing the lower tail of interest.

## Reading

The census sharpens the successor-paper framing from the discussion: over five decades of n, the multiplicity distribution of the corner strip behaves as the block-random model predicts in every statistic measured — mean pinned to H_N, PNT rate at R = 1, P₂ stratum at the log log law, no square-centered anomaly above noise — while the one systematic deviation is the universal 2e^(−γ) factor separating the independence model from the truth at exactly the stratum a proof must control. The empirical content of the Jₙ hypothesis is thus isolated: it is not that J_n is atypical (it is not, to n^(−1/2) resolution), but that typicality at the R = 1 stratum is exactly what no unconditional method certifies. The numbers put the gap where the theory says it is, and nowhere else.

## Files

`rn_census.py` — preregistered script (predictions in header, run unmodified). `summary.csv` — anchor-level statistics. `histograms.zip` — full R-histograms, both windows, all anchors. `chart_distribution.png` — multiplicity distribution at n = 10⁷, both windows, log-log. `chart_ratios.png` — P2 and P3 convergence with the 2e^(−γ) line. `chart_tv.png` — P4 total-variation decay against the n^(−1/2) guide.

Runtime: 85 s total, single-threaded, peak memory ~160 MB (n = 10⁷).
