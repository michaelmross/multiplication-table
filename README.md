# Per-gap table density at astronomic heights: run protocol

## What this campaign is

Three data sources bracketing one quantity — the density ρ(N, c) of integers near (cN)² with a divisor within factor 1/c of their square root, equivalently the distinct-entry density of the N × N multiplication table in the square-gap at k = cN:

1. **Exact census** (done): `gap_summary.csv`, N ≤ 10⁸, exact divisor sieve with brute-force cross-check.
2. **Tier 1, exact integers** (yours to run): `t1_exact_mc.py`, Kalai random factored integers, c = 1/2, N from 10⁸ (validation against the exact value 0.25276) up to 2^1024. Anchors the model to true arithmetic.
3. **Tier 2, continuum model** (done to 2^(10^11); extendable): `t2_model_mc.py`, GEM(1) stick-breaking factorizations, three window widths c ∈ {0.5, 0.75, 0.9}, heights unbounded because only L = ln x enters. Existing data: `mc_tier2.csv`, 41 rows, t = log log x from 2.8 to 25.7.

Preregistered predictions P1–P5 and their evaluation criteria are frozen in the script headers and implemented verbatim in `analyze_mc.py`. Provisional verdicts on the existing Tier-2 data: P2 fails at narrow windows (model deficit up to 23% at c = 0.9); P3 shape emergence is visible but the strict all-intervals criterion narrowly fails (two intervals at t ≈ 9–13 deviate by 0.05–0.07); P4 universality arrives near t ≈ 16 rather than 8; P5 pooled d̂ ≈ 0.05 sits below the 0.086 ± 0.03 band. Your runs can move P1 (unrun), the seam offsets, and — with more samples at the giants — P5.

## What to run

Dependencies: Python 3.10+, numpy, gmpy2 (`pip install gmpy2`).

**Tier 1** (the piece only exact integers can give):

    python3 t1_exact_mc.py --min-hits 400 --procs <cores> --hours 12

Runs N = 10⁸, 2^64, 2^128, 2^256, 2^512, 2^1024 in order, appending to `mc_tier1.csv`. Each height starts with a 30 s single-core warm-up that prints attempts/sample and an ETA before committing — kill and rerun with a smaller `--hours` or stop after 2^256 if the top heights are dearer than expected. Cost is dominated by gmpy2 primality tests on chain elements; expect roughly seconds per accepted sample per core at the 2^512–2^1024 heights. `--min-hits 2500` upgrades relative SE from ~5% to ~2% at proportional cost; the P1 validation and the 2^64/2^256 seams are the priority, the 2^1024 seam is the luxury.

**Tier 2 extension** (optional, where cycles buy the most inference):

    python3 t2_model_mc.py --target-hits 2500 --tcap 7200 --procs <cores>

Re-runs the full height grid at ~2% relative SE including 2^(10^12) (t = 27.96), appending to `mc_tier2.csv`. The P5 δ-recovery fit is currently variance-limited exactly at the giant heights (ρ ~ 0.002–0.008 there), so this is the run that decides whether d̂ climbs from ~0.05 toward 0.086 or the shortfall is real — either answer is interesting, since the GEM model is literally the uniform-order-statistics mechanism Ford's constants come from, and a persistent shortfall at t ≈ 28 would mean even the *model's* asymptopia lies deeper still.

**Analysis** (run in the directory containing whichever CSVs exist):

    python3 analyze_mc.py

Prints PASS/FAIL per criterion plus the seam offsets, and writes `mc_chart_slopes.png` (the payoff chart: local decay rate against Ford's −δ − (3/2)/t across 23+ e-folds of t, exact census + model + your Tier-1 points) and `mc_chart_density.png`.

## Interrupt and resume

Both runners are safe to kill and restart against the same `--out` CSV. Tier 1 commits work in rounds (default 20 min, `--round-mins`): each round appends its own row under the same tag, so a kill forfeits at most one round; on restart, banked hits are summed per tag and any height already at `--min-hits` is skipped (`--redo` overrides). The warm-up's own samples are banked too. Tier 2 commits per (height, c) cell — cells are minutes-scale, so a kill loses at most one cell, and completed cells are skipped on restart. Seeds advance per round and per cell, so resuming never replays an RNG stream. `analyze_mc.py` pools duplicate-tag rows by sample count before computing anything, so round-fragmented data analyzes identically to a single run. The one thing resume does *not* protect: appending rows from runs with different physical parameters (a changed height grid, a changed c set) into the same CSV — keep one CSV per campaign configuration. Relatedly, if you start a run in a **fresh directory** (no CSV present) and later concatenate its rows with an earlier run's, use a different `--seed` for the second run: with the default seed, both sessions begin from the same RNG streams, so their early samples are identical rather than independent, and pooled standard errors become optimistic. Within a single directory this cannot happen — round and cell seeds now advance with the number of rows already banked.

## Honesty notes

Seeds are baked in but changing them is fine — the preregistration binds criteria, not RNG streams. The Tier-2 window is centered at M/2 (sum of kept sticks), correcting an L/2 mis-centering caught during profiling before any measurement data was taken; all rows in `mc_tier2.csv` postdate the fix. The model's known deficits (no prime multiplicities, no small-prime discreteness) are measured at the seams, not assumed away — that is what Tier 1 is for. If a Tier-1 height finishes with n < 200 kept samples, treat its seam offset as indicative only.

## Files

`t1_exact_mc.py`, `t2_model_mc.py`, `analyze_mc.py`, `mc_tier2.csv` (41 completed rows), plus the earlier deliverables `gap_census.py`, `gap_summary.csv`, `rn_census.py`, `summary.csv`.

## Campaign A-prime (constant floor vs 1/t tail)

`campaign_a2.py run` tops up every cell from 2^10^8 through the new arm at
2^10^16, 2^10^17, 2^10^18 (t = 37.2, 39.5, 41.8) to frozen hit targets,
reading banked hits from mc_tier2.csv and running only shortfalls — safe to
kill and relaunch indefinitely. Frozen constants: C_frozen = 0.0126,
b_frozen = 0.359. `campaign_a2.py evaluate` renders D1 (window precisions),
D2 (deep-arm z against both frozen hypotheses), and D3 (two-parameter fit
deficit = C_inf + b/t; a 3-sigma nonzero C_inf means the model's asymptotic
exponent is not Ford's delta). Cost at default targets: ~130-160 CPU-hours
total at 2-3 ms/sample; --scale multiplies all targets. The 2^10^18 ceiling
is set by a sparsity argument for longdouble validity (logged in the script
header); beyond it lies double-double territory.

## Campaign A3 (exact-window deep arm)

`campaign_a3.py validate` first: a paired exact-vs-longdouble comparison on
identical stick streams. The clean-height gate must pass; the 60k-pair deep
runs are the ARBITER between the two live explanations of the deep-arm
anomaly — nonzero one-sided flips confirm the longdouble blur artifact and
its sign, while ~zero flips vindicate longdouble and mean the anomaly may be
real physics. Either way the exact engine is correct, so proceed to
`campaign_a3.py run` (writes mc_tier3.csv only; the voided tier-2 rows for
2^10^17/2^10^18 are never counted) and finish with `campaign_a3.py
evaluate`, which renders the frozen D1/D2/D3 on the merged clean set. The
exact representation is valid to L < 2^62; expect ~5-9 ms/sample and
~20-30 pool-hours for the four cells at default targets.
