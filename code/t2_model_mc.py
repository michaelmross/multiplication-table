#!/usr/bin/env python3
"""
Tier-2: continuum (GEM) Monte Carlo, consolidated final version, standalone.

Model. For a uniform integer of magnitude ~ x, the log-prime multiset
conditioned on its total is L * GEM(1) stick-breaking, L = ln x (Billingsley/
PD(1)).  Sticks below ln 2 are discarded (no primes below 2); the model
integer's log-size is M = sum of kept sticks and the divisor window is
[M/2 - h, M/2 + h], h = ln(1/c).  Known model limitations (measured at the
seams, not assumed away): no prime multiplicities, no small-prime
discreteness -- the measured deficit vs exact truth is ~8% (c=0.5) to ~23%
(c=0.9) at accessible heights, shrinking with height at c=0.5.

Subset test: pruned meet-in-the-middle; DFS branch-and-bound fallback for
stick counts > 36 (fast when hits are plentiful, which is the large-K regime).

Usage:
    python3 t2_model_mc.py [--target-hits H] [--tcap S] [--procs P] [--out CSV]
Defaults H=600 (rel SE ~4%); for the delta-recovery fit (P5) at the giant
heights use H=2500 and tcap 3600+.  Heights: overlap {10^4,10^6,10^8},
2^64 .. 2^2^20, then 2^10^8, 2^10^10, 2^10^11, 2^10^12.
c in {0.5, 0.75, 0.9}.  Appends: tag, c, L, t, rho, se, n, seconds.
"""
import argparse, csv, math, os, random, time
import multiprocessing as mp
import numpy as np

LN2 = math.log(2.0)

def hit_mim(vals, w1, w2):
    A = np.zeros(1); B = np.zeros(1)
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
    suf = np.concatenate((np.cumsum(vals[::-1])[::-1], [0.0]))
    K = len(vals)
    stack = [(0, 0.0)]
    while stack:
        i, s = stack.pop()
        if w1 <= s <= w2:
            return True
        if i == K or s > w2 or s + suf[i] < w1:
            continue
        stack.append((i + 1, s))
        stack.append((i + 1, s + vals[i]))
    return False

def hit(vals, w1, w2):
    vals = [v for v in vals if v <= w2]
    if not vals or sum(vals) < w1:
        return False
    vals.sort(reverse=True)
    for v in vals:
        if w1 <= v <= w2:
            return True
    if len(vals) > 36:
        return hit_dfs(np.array(vals), w1, w2)
    return hit_mim(vals, w1, w2)

def one(L, h, rng):
    rem = L; s = []
    while rem >= LN2:
        w = rem * rng.random()
        if w >= LN2:
            s.append(w)
        rem -= w
    M = sum(s)
    return hit(s, M/2 - h, M/2 + h)

def worker(args):
    L, h, n, seed, tcap = args
    rng = random.Random(seed)
    t0 = time.time(); hits = tot = 0
    while tot < n and time.time() - t0 < tcap:
        for _ in range(min(50, n - tot)):
            hits += one(L, h, rng); tot += 1
    return hits, tot

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-hits", type=int, default=600)
    ap.add_argument("--tcap", type=float, default=600.0)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--out", default="mc_tier2.csv")
    ap.add_argument("--seed", type=int, default=20260823)
    ap.add_argument("--redo", action="store_true",
                    help="rerun cells even if rows already banked")
    a = ap.parse_args()
    heights = [
        ("N=10^4", math.log(10**4)), ("N=10^6", math.log(10**6)),
        ("N=10^8", math.log(10**8)),
        ("2^64", 64*LN2), ("2^256", 256*LN2), ("2^1024", 1024*LN2),
        ("2^4096", 4096*LN2), ("2^16384", 16384*LN2), ("2^65536", 65536*LN2),
        ("2^2^18", 2**18*LN2), ("2^2^20", 2**20*LN2),
        ("2^10^8", 1e8*LN2), ("2^10^10", 1e10*LN2),
        ("2^10^11", 1e11*LN2), ("2^10^12", 1e12*LN2),
    ]
    done = {}
    if os.path.exists(a.out):
        with open(a.out) as g:
            for r in csv.DictReader(g):
                key = (r["tag"], float(r["c"]))
                done[key] = done.get(key, 0) + int(r["n"])
    new = not os.path.exists(a.out)
    f = open(a.out, "a", newline=""); w = csv.writer(f)
    if new:
        w.writerow(["tag","c","L","t","rho","se","n","seconds"])
    si = 0
    for tag, lnN in heights:
        for c in (0.5, 0.75, 0.9):
            si += 1
            if done.get((tag, c), 0) > 0 and not a.redo:
                print(f"{tag:>8} c={c:.2f}  {done[(tag, c)]} samples banked, "
                      f"skipping (--redo to override)", flush=True)
                continue
            L = 2*lnN + 2*math.log(c); h = math.log(1/c)
            t0 = time.time()
            # pilot on one core
            prior = done.get((tag, c), 0)
            rng = random.Random(a.seed + si + 7919 * prior)
            ph = sum(one(L, h, rng) for _ in range(2000))
            rho0 = max(ph/2000, 1e-9)
            n_target = min(2_000_000, max(2000, int(a.target_hits / rho0)))
            per = -(-n_target // a.procs)
            args = [(L, h, per, a.seed + si*10**4 + 7919*prior + j + 1,
                     a.tcap) for j in range(a.procs)]
            with mp.Pool(a.procs) as pool:
                res = pool.map(worker, args)
            hits = ph + sum(r[0] for r in res)
            n = 2000 + sum(r[1] for r in res)
            rho = hits/n; se = math.sqrt(rho*(1-rho)/n)
            w.writerow([tag, c, L, math.log(L), rho, se, n,
                        round(time.time()-t0, 1)]); f.flush()
            print(f"{tag:>8} c={c:.2f}  t={math.log(L):6.2f}  "
                  f"rho={rho:.5f}+-{se:.5f}  n={n}", flush=True)
    f.close()

if __name__ == "__main__":
    main()
