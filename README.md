# The Multiplication Table Near Perfect Squares

Computational examination of the Erdős–Tenenbaum–Ford constant
δ = 1 − (1 + log log 2)/log 2 = 0.086071…, localized to square-gaps of the
N × N multiplication table. This repository contains the complete code,
data, preregistration ledger, and paper source for the note [*The
Multiplication Table Near Perfect Squares: Concentration, Typicality, and
the Slow Approach to Ford's Exponent*](https://zenodo.org/records/22286360) (Ross, M. M., 2026).

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.1354181135-blue.svg)](https://doi.org/10.5281/zenodo.22287094)

## Findings

The quantity ρ(N, c) is the density of distinct N × N-table entries in the
square-gap ((cN)², (cN+1)²], equivalently the density of integers near
(cN)² with a divisor within factor 1/c of their square root. Measured by
exact sieve to N = 10⁸, exact-integer Monte Carlo to N = 2^1024, and a
seam-anchored continuum model of factorizations to model heights of
N = 2^(10^18):

- **Typicality.** Square-centered windows are statistically
  indistinguishable from generic windows of the same magnitude
  (powered ensemble null: z = −0.02, +0.90, −0.26 at N = 10⁶, 10⁷, 10⁸).
- **Concentration.** The gap-to-gap coefficient of variation fits N^β with
  β = −0.440 ± 0.031. The unproven localization of Ford's theorem to
  √x-scale windows is empirically a concentration statement with a
  square-root exponent.
- **Shape and exponent.** The local decay attains Ford's shape
  −δ − (3/2)/t in t = log log x only from t ≈ 13–16. The exponent itself
  is never attained at any computed height. The deficit follows a 1/t law
  with coefficient b = 0.333 ± 0.032 in the mature regime, a constant
  floor being rejected at Δχ² = 56.8, behind a transient dying near
  t ≈ 25. At t = 39.5, the model height of quintillion-bit integers, the
  effective exponent still reads 0.077.
- **Mechanism.** An exactly computable ballot skeleton (a killed
  negative-drift walk, `code/ballot_skeleton.py`, milliseconds to run)
  reproduces δ, the t^(−3/2), and the 1/t coefficient to 1–2 standard
  errors. The subleading law is a property of the order-statistics
  mechanism itself. The skeleton predicts b declines toward ~0.22 by
  t ≈ 50, a falsifiable target for future campaigns.

No quintillion-bit integer was ever represented. The double logarithm
compresses that height to t ≈ 40, and a model sample is about forty stick
values. The largest genuine integers computed are the 2046-bit numbers of
the Kalai tier. See the note, Section 2.3.

## Repository map

    README.md               this file
    paper/                  note.tex, note.pdf, and the figures it embeds
    code/                   all campaign scripts (inventory below)
    results/                all data tiers and campaign reports
    figures/                all charts, including superseded ones

Reports: `results/FINAL_report.md` is the campaign chronicle and full
preregistration ledger, including every refuted prediction and both
instrument errata. `results/RN_census_report.md` and
`results/gap_census_report.md` are the exact-height censuses.

## Reproduction guide

Requirements: Python ≥ 3.10, numpy, gmpy2, matplotlib. Campaigns A and A′
require genuinely extended-precision `numpy.longdouble` and abort
otherwise, so run under Linux or WSL, not Windows-native Python. Scripts
import from each other (`campaign_a2.py` and `campaign_a3.py` import
`campaign_a.py`), so keep `code/` together.

Pipeline, in dependency order, with rough costs on a modern multicore
machine:

1. `rn_census.py` — strip-multiplicity census over J_n (minutes).
2. `gap_census.py` — exact per-gap density census, brute-force verified
   (under a minute).
3. `campaign_b.py ensemble` — powered typicality null and concentration
   exponent (about an hour).
4. `t1_exact_mc.py` — exact-integer ladder to 2^1024 via Kalai sampling
   (hours to days by `--min-hits`; resumable in rounds).
5. `t2_model_mc.py` — continuum model across the height grid (hours).
6. `campaign_a.py` — deep-t extension to 2^(10^15) in extended precision
   (overnight).
7. `campaign_a2.py run` — precision boost plus the 2^(10^16..18) arm
   (~130–160 CPU-hours). The 2^(10^17) and 2^(10^18) cells it produces
   are **voided** (see below) and retained for the ledger.
8. `campaign_a3.py validate`, then `run`, then `evaluate` — the paired
   exact-vs-extended arbiter, the exact-engine rerun of the deep arm
   (~20–30 pool-hours), and the frozen final verdicts.
9. `ballot_skeleton.py` — the mechanism companion (milliseconds).
10. `analyze_mc.py` — mid-campaign evaluation and charts.

Every quantitative prediction was frozen in the script headers before its
data existed. Evaluation criteria are implemented verbatim in
`analyze_mc.py`, `campaign_b.py evaluate`, `campaign_a2.py evaluate`, and
`campaign_a3.py evaluate`. Data files append-only with resume support.
Seeds advance with banked rows, so interrupted runs never replay an RNG
stream. If you start a fresh directory and later merge CSVs, use a
different `--seed` (see the note's Appendix A and the script headers).

## Data integrity notes

`results/mc_tier2.csv` is the final file, including the Campaign A′ boost
rows and the voided extended-precision rows for 2^(10^17) and 2^(10^18).
The voided rows are deliberately retained. The evaluation code excludes
them by tag, and the paired arbiter that voided them (seven classification
flips in 240,000 pairs, all one-sided, matching the observed anomaly to
Poisson precision) is documented in the note's Appendix A. The exact-engine
replacements live in `results/mc_tier3.csv`.

## Citation and license

Code: MIT. Note and data: CC-BY 4.0. If you use this work, cite the
Zenodo record (DOI:10.5281/zenodo.22286360).
