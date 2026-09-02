#!/usr/bin/env python3
"""
Campaign A: the exponent deficit at t = 28 .. 35  (N up to 2^(10^15)).

Question. The completed campaign measured the effective exponent
d_hat(t) = -(d ln rho/dt + 3/2 t^{-1}) climbing toward delta = 0.086071
but still short by 0.0122 +/- 0.0055 pooled over t in [20, 28.5].  Is the
deficit a genuine subleading correction, and does it decay like 1/t?

FROZEN BASELINE (computed from the completed mc_tier2.csv BEFORE this run;
do not refit):
    pooled deficit on t in [20, 28.5]:  +0.0122 +/- 0.0055
    1/t fit through origin:             b_hat = 0.278 +/- 0.123
    predicted deficit at t = 31.5:      +0.0088 under H_{1/t};  0 under H_0.

PREREGISTERED EVALUATION (implemented verbatim in analyze_deep.py):
  A1. New heights 2^10^13, 2^10^14, 2^10^15 at c in {0.5, 0.75} with
      >= 4000 hits per cell yield a pooled deficit on t in [28, 35.5] with
      SE <= 0.005.  The pooled value is reported with its SE regardless.
  A2. Three-way decision at the frozen prediction:
      z0 = deficit/SE (against H_0: already asymptotic, deficit 0),
      z1 = (deficit - 0.278/t_bar)/SE (against H_{1/t} with frozen b_hat).
      Outcomes: excludes H_0 (|z0| > 2, |z1| <= 2); excludes H_{1/t}
      (|z1| > 2, |z0| <= 2); excludes both; distinguishes neither.
      A joint weighted fit of b over ALL windows t >= 20 (old + new) is
      also reported.
  Design choices made in advance: c = 0.9 is dropped (precision-marginal
  window, sample-hungriest, universality already established); height
  ceiling 2^(10^15) is set by extended-precision headroom (see below).

PRECISION.  All stick generation, subset sums, and window tests use
numpy longdouble (80-bit extended on x86-64 Linux: eps ~ 1.1e-19, ulp at
L = 1.4e16 is ~1.5e-3, rms accumulated error over K ~ 37 sticks ~ 0.01,
against window half-widths h >= 0.288).  The script ABORTS if longdouble
is not genuinely extended (Windows-native numpy aliases it to float64 --
run under WSL).

Resume: identical CSV/row protocol to t2_model_mc.py; appends to
mc_tier2.csv; completed cells are skipped; seeds advance with banked rows.

Usage:  python3 campaign_a.py [--target-hits 4000] [--tcap 999999]
                              [--procs P] [--out mc_tier2.csv]
Cost guidance: mean stick count K ~ 35-37 at the top heights, so
meet-in-the-middle sides reach ~2^18; expect tens of ms per sample and
plan for overnight-scale accumulation per top cell.  Rounds are cheap to
interrupt; relaunches resume.
"""
import argparse, csv, math, os, random, time
import multiprocessing as mp
import numpy as np

LD = np.longdouble
LN2 = LD(math.log(2.0))

def precision_guard():
    eps = np.finfo(np.longdouble).eps
    if eps > 1e-18:
        raise SystemExit(
            f"numpy longdouble eps = {eps} (not extended precision). "
            "On Windows-native Python longdouble aliases float64 -- "
            "run this under WSL.")

def hit_mim(vals, w1, w2):
    A = np.zeros(1, dtype=LD); B = np.zeros(1, dtype=LD)
    for i, v in enumerate(vals):
        if i % 2 == 0:
            A = np.concatenate((A, A + v)); A = A[A <= w2]
        else:
            B = np.concatenate((B, B + v)); B = B[B <= w2]
    B.sort()
    lo = np.searchsorted(B, w1 - A, side="left")
    hi = np.searchsorted(B, w2 - A, side="right")
    return bool(np.any(hi > lo))

def hit_dfs(vals, w1, w2):
    arr = np.array(vals, dtype=LD)
    suf = np.concatenate((np.cumsum(arr[::-1])[::-1],
                          np.zeros(1, dtype=LD)))
    K = len(arr)
    stack = [(0, LD(0.0))]
    while stack:
        i, s = stack.pop()
        if w1 <= s <= w2:
            return True
        if i == K or s > w2 or s + suf[i] < w1:
            continue
        stack.append((i + 1, s))
        stack.append((i + 1, s + arr[i]))
    return False

def hit(vals, w1, w2):
    vals = [v for v in vals if v <= w2]
    if not vals:
        return False
    tot = LD(0.0)
    for v in vals:
        tot += v
    if tot < w1:
        return False
    vals.sort(reverse=True)
    for v in vals:
        if w1 <= v <= w2:
            return True
    if len(vals) > 36:
        return hit_dfs(vals, w1, w2)
    return hit_mim(vals, w1, w2)

def one(L, h, rng):
    rem = L
    s = []
    while rem >= LN2:
        w = rem * LD(rng.random())
        if w >= LN2:
            s.append(w)
        rem -= w
    M = LD(0.0)
    for v in s:
        M += v
    return hit(s, M / 2 - h, M / 2 + h)

def worker(args):
    L, h, n, seed, tcap = args
    rng = random.Random(seed)
    t0 = time.time(); hits = tot = 0
    while tot < n and time.time() - t0 < tcap:
        for _ in range(min(20, n - tot)):
            hits += one(L, h, rng); tot += 1
    return hits, tot

def main():
    precision_guard()
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-hits", type=int, default=4000)
    ap.add_argument("--tcap", type=float, default=999999.0)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--out", default="mc_tier2.csv")
    ap.add_argument("--seed", type=int, default=20260824)
    a = ap.parse_args()
    heights = [("2^10^13", 1e13 * math.log(2.0)),
               ("2^10^14", 1e14 * math.log(2.0)),
               ("2^10^15", 1e15 * math.log(2.0))]
    done = {}
    if os.path.exists(a.out):
        with open(a.out) as g:
            for r in csv.DictReader(g):
                key = (r["tag"], float(r["c"]))
                done[key] = done.get(key, 0) + 1
    new = not os.path.exists(a.out)
    f = open(a.out, "a", newline=""); w = csv.writer(f)
    if new:
        w.writerow(["tag","c","L","t","rho","se","n","seconds"])
    si = 0
    for tag, lnN in heights:
        for c in (0.5, 0.75):
            si += 1
            prior = done.get((tag, c), 0)
            # one CSV row per launch-session per cell; resume = relaunch
            Lf = 2 * lnN + 2 * math.log(c)
            L = LD(Lf); h = LD(math.log(1.0 / c))
            t0 = time.time()
            rng = random.Random(a.seed + si + 7919 * prior)
            pilot = 1000
            ph = sum(one(L, h, rng) for _ in range(pilot))
            rho0 = max(ph / pilot, 1e-9)
            n_target = min(5_000_000, max(pilot,
                           int(a.target_hits / rho0)))
            per = -(-n_target // a.procs)
            args = [(L, h, per,
                     a.seed + si * 10**4 + 7919 * prior + j + 1, a.tcap)
                    for j in range(a.procs)]
            with mp.Pool(a.procs) as pool:
                res = pool.map(worker, args)
            hits = ph + sum(r[0] for r in res)
            n = pilot + sum(r[1] for r in res)
            rho = hits / n; se = math.sqrt(rho * (1 - rho) / n)
            w.writerow([tag, c, Lf, math.log(Lf), rho, se, n,
                        round(time.time() - t0, 1)]); f.flush()
            print(f"{tag} c={c}  t={math.log(Lf):.2f}  "
                  f"rho={rho:.5f}+-{se:.5f}  n={n}  "
                  f"[{(time.time()-t0)/3600:.2f} h]", flush=True)
    f.close()

if __name__ == "__main__":
    main()
