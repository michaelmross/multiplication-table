#!/usr/bin/env python3
"""
Campaign B: ensemble exact census (retires the single-window caveat), plus
the preregistered A1/A2 evaluation for Campaign A -- run with subcommands:

    python3 campaign_b.py ensemble [--procs P]     # exact-sieve ensembles
    python3 campaign_b.py evaluate                 # A1/A2 on mc_tier2.csv

ENSEMBLE (preregistered):
  For N in {10^6, 10^7, 10^8}:
  B1. Multiplicity histograms R_N over J_n and 32 disjoint shifted windows
      [4n^2 + 5jn, 4n^2 + 5jn + 2n], j = 1..32, n = N/2.  Criterion: the
      TV distance of J_n's histogram from the ensemble mean lies within
      z <= 2 of the shifted windows' own TV-to-mean distribution.
      (Upgrades the earlier one-window diagnostic to a powered test of
      square-centered typicality.)
  B2. Distinct N x N-table entries D over 50 consecutive square-gaps at
      k = N/2.  Criterion: fitting CV(D) ~ N^beta across the three anchors
      (plus the earlier 10-gap values as a cross-check) gives
      beta = -0.5 +/- 0.1.
  Output: ensemble_summary.csv, per-anchor histogram files.

EVALUATE (implements Campaign A's frozen criteria verbatim):
  Baseline frozen from the pre-Campaign-A mc_tier2.csv:
      pooled deficit on t in [20, 28.5]:  +0.0122 +/- 0.0055
      1/t coefficient:                    b_hat = 0.278 +/- 0.123
  A1. Pooled deficit delta - d_hat over interval midpoints t in [28, 35.5]
      (c in {0.5, 0.75}), reported with SE; target SE <= 0.005.
  A2. z0 = deficit/SE against H_0 (deficit 0);
      z1 = (deficit - 0.278/t_bar)/SE against H_{1/t}.
      Verdict: excludes H_0 / excludes H_{1/t} / both / neither.
      Also reports the joint weighted 1/t fit over all t >= 20.
"""
import argparse, csv, math, os, time
import multiprocessing as mp
import numpy as np

DELTA = 1.0 - (1.0 + math.log(math.log(2.0))) / math.log(2.0)
B_FROZEN = 0.278          # 1/t coefficient, frozen pre-Campaign-A
DEF_FROZEN = 0.0122       # pooled deficit t in [20, 28.5], frozen

# ---------------- ensemble ----------------
def strip_hist(args):
    L, W, N = args
    R = np.ones(W, dtype=np.int32)
    for a in range(2, N + 1):
        R[(-L) % a::a] += 1
    return np.bincount(R)

def gap_D(args):
    N, k = args
    lo = k * k + 1; hi = (k + 1) * (k + 1)
    hit = np.zeros(hi - lo + 1, dtype=bool)
    a = max(1, -(-lo // N))
    while a <= N:
        m = ((lo + a - 1) // a) * a
        cap = min(hi, a * N)
        while m <= cap:
            hit[m - lo] = True
            m += a
        a += 1
    return int(hit.sum()), hi - lo + 1

def gap_D_fast(args):
    # vectorized version of gap_D (same result, numpy chunks)
    N, k = args
    lo = k * k + 1; hi = (k + 1) * (k + 1)
    hit = np.zeros(hi - lo + 1, dtype=bool)
    a_lo = max(1, -(-lo // N))
    while a_lo <= N:
        a_hi = min(N, a_lo + 5_000_000 - 1)
        a = np.arange(a_lo, a_hi + 1, dtype=np.int64)
        m = ((lo + a - 1) // a) * a
        cap = np.minimum(np.int64(hi), a * np.int64(N))
        while True:
            ok = m <= cap
            if not ok.any():
                break
            hit[m[ok] - lo] = True
            m = m + a
        a_lo = a_hi + 1
    return int(hit.sum()), hi - lo + 1

def tv(h1, h2):
    k = max(len(h1), len(h2))
    a = np.zeros(k); b = np.zeros(k)
    a[:len(h1)] = h1; b[:len(h2)] = h2
    return 0.5 * np.abs(a / a.sum() - b / b.sum()).sum()

def run_ensemble(procs):
    out = open("ensemble_summary.csv", "w", newline="")
    w = csv.writer(out)
    w.writerow(["N", "tv_Jn_z", "tv_Jn", "tv_shift_mean", "tv_shift_sd",
                "cv_D_50gaps", "seconds"])
    for N in (10**6, 10**7, 10**8):
        n = N // 2
        t0 = time.time()
        jobs = [(4*n*n - n, 2*n + 1, N)] + \
               [(4*n*n + 5*j*n, 2*n + 1, N) for j in range(1, 33)]
        with mp.Pool(procs) as pool:
            hists = pool.map(strip_hist, jobs)
        kmax = max(len(h) for h in hists)
        padded = np.array([np.pad(h, (0, kmax - len(h))) for h in hists],
                          dtype=float)
        mean_shift = padded[1:].mean(axis=0)
        tv_j = tv(padded[0], mean_shift)
        tv_s = np.array([tv(padded[i], mean_shift)
                         for i in range(1, 33)])
        z = (tv_j - tv_s.mean()) / tv_s.std(ddof=1)
        with mp.Pool(procs) as pool:
            Ds = pool.map(gap_D_fast, [(N, n + j) for j in range(50)])
        D = np.array([d for d, _ in Ds], dtype=float)
        cv = D.std(ddof=1) / D.mean()
        w.writerow([N, round(z, 3), tv_j, tv_s.mean(), tv_s.std(ddof=1),
                    cv, round(time.time() - t0, 1)])
        out.flush()
        np.savetxt(f"ensemble_hist_N{N}.csv", padded, delimiter=",")
        print(f"N={N}: B1 z(TV of J_n vs ensemble) = {z:+.2f}  "
              f"B2 CV(D, 50 gaps) = {cv:.5f}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    out.close()
    print("B1 verdict per anchor: PASS iff |z| <= 2.")
    print("B2: fit CV ~ N^beta across anchors; PASS iff "
          "beta in [-0.6, -0.4].")

# ---------------- evaluate ----------------
def run_evaluate():
    rows = [dict(r) for r in csv.DictReader(open("mc_tier2.csv"))]
    acc = {}
    for r in rows:
        key = (r["tag"], float(r["c"]))
        n = int(r["n"]); h = float(r["rho"]) * n
        if key in acc:
            acc[key]["n"] += n; acc[key]["h"] += h
        else:
            acc[key] = {"t": float(r["t"]), "n": n, "h": h}
    pts = []
    for c in (0.5, 0.75):
        sub = sorted((v | {"c": c} for (tag, cc), v in acc.items()
                      if cc == c), key=lambda v: v["t"])
        t = np.array([v["t"] for v in sub])
        rho = np.array([v["h"] / v["n"] for v in sub])
        sig = np.sqrt((1 - rho) / (rho * np.array([v["n"] for v in sub])))
        tm = (t[1:] + t[:-1]) / 2
        sl = np.diff(np.log(rho)) / np.diff(t)
        se = np.sqrt(sig[1:]**2 + sig[:-1]**2) / np.diff(t)
        for a, b, s in zip(tm, sl, se):
            pts.append((a, DELTA - (-(b + 1.5 / a)), s))
    new = [(a, d, s) for a, d, s in pts if 28 <= a <= 35.5]
    if not new:
        print("no intervals with midpoint in [28, 35.5] yet -- "
              "run campaign_a.py first")
        return
    t_, d_, s_ = map(np.array, zip(*new))
    wgt = 1 / s_**2
    deficit = float(np.sum(d_ * wgt) / np.sum(wgt))
    se_p = 1 / math.sqrt(np.sum(wgt))
    tbar = float(np.sum(t_ * wgt) / np.sum(wgt))
    z0 = deficit / se_p
    z1 = (deficit - B_FROZEN / tbar) / se_p
    print(f"A1: pooled deficit on t in [28, 35.5] = {deficit:+.4f} "
          f"+- {se_p:.4f} (target SE <= 0.005: "
          f"{'MET' if se_p <= 0.005 else 'NOT MET'})")
    print(f"A2: t_bar = {tbar:.1f}; H_1/t predicts "
          f"{B_FROZEN/tbar:+.4f}; z0 = {z0:+.2f} (vs deficit 0), "
          f"z1 = {z1:+.2f} (vs frozen 1/t)")
    if abs(z0) > 2 and abs(z1) <= 2:
        print("A2 verdict: excludes H_0; consistent with 1/t decay.")
    elif abs(z1) > 2 and abs(z0) <= 2:
        print("A2 verdict: excludes H_1/t; consistent with asymptotic.")
    elif abs(z0) > 2 and abs(z1) > 2:
        print("A2 verdict: excludes BOTH frozen hypotheses.")
    else:
        print("A2 verdict: does not distinguish the hypotheses.")
    old = [(a, d, s) for a, d, s in pts if a >= 20]
    t_, d_, s_ = map(np.array, zip(*old))
    b = float(np.sum(d_ * (1/t_) / s_**2) / np.sum((1/t_)**2 / s_**2))
    bse = 1 / math.sqrt(float(np.sum((1/t_)**2 / s_**2)))
    print(f"joint 1/t fit over all t >= 20 (old + new): "
          f"b = {b:+.3f} +- {bse:.3f} (frozen baseline 0.278)")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["ensemble", "evaluate"])
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    a = ap.parse_args()
    if a.mode == "ensemble":
        run_ensemble(a.procs)
    else:
        run_evaluate()
