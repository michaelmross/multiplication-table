#!/usr/bin/env python3
"""
Tier-1: exact-integer Monte Carlo for the per-gap table density, standalone.

Quantity. rho(N, c=1/2) = density of integers m ~ x = (N/2)^2 having a
divisor within factor 2 of sqrt(m).  Matches the square-gap distinct-entry
density of the exact census (k = N/2) and the Tier-2 continuum model at the
same heights (L = ln x), giving the two seams:
    exact census (N <= 10^8)  <->  Tier-1 (10^8 .. 2^1024)  <->  Tier-2 (all)

Method: Kalai's algorithm for uniform random factored integers <= x
(chain x >= s_1 >= s_2 >= ... -> 1 with s_{i+1} uniform in [1, s_i];
m = product of prime s_i with repetition; accept iff m <= x w.p. m/x, both
exact in integer arithmetic).  Keep samples with m > x/e (63%; density drift
across one e-fold of ln m is O(1/ln x) in t and negligible).  Divisor test:
subset of the factor multiset's logs summing into
[ln(m)/2 - ln 2, ln(m)/2 + ln 2].  Early abort when the partial product
exceeds x.

Preregistered evaluation (see README_mc.md):
  P1: at N = 10^8 this reproduces the exact census value rho = 0.25276
      within max(2 SE, 1%).
  Seam: at N in {2^64, 2^256, 2^1024} the Tier-2 model offset
      (model/exact - 1) is measured; P2 bound: |offset| <= 20%.

INTERRUPT/RESUME: work is committed in rounds (default 20 min); each round
appends its own CSV row (same tag), so a kill loses at most one round.  On
restart, rows already in --out are read, banked hits are summed per tag, and
heights with banked hits >= --min-hits are skipped (override with --redo).
analyze_mc.py pools duplicate-tag rows by sample count, so partial rounds
combine correctly.  Seeds advance per round; resuming never reuses a stream.

Usage:
    python3 t1_exact_mc.py [--min-hits H] [--procs P] [--hours T]
                           [--round-mins R] [--out CSV] [--redo]
Defaults: H=400 (rel SE ~5%; use 2500 for ~2%), P=all cores, T=12.0 per
height, R=20.  Heights: N = 10^8, 2^64, 2^128, 2^256, 2^512, 2^1024.
A warm-up at each fresh height prints attempts/sample and an ETA (and banks
its own samples as a row).  Requires: gmpy2.
"""
import argparse, csv, math, os, random, time
import multiprocessing as mp
import numpy as np
from gmpy2 import is_prime as _gmp_is_prime

LN2 = math.log(2.0)
_SS_LIM = 10_000_000
_sieve = None

def _init_sieve():
    global _sieve
    if _sieve is None:
        s = np.ones(_SS_LIM, dtype=bool); s[:2] = False
        for p in range(2, 3163):
            if s[p]:
                s[p*p::p] = False
        _sieve = s

def is_prime_int(s):
    if s < _SS_LIM:
        return bool(_sieve[s])
    return bool(_gmp_is_prime(s))

def subset_hit(logs, w1, w2):
    logs = [v for v in logs if v <= w2]
    if not logs or sum(logs) < w1:
        return False
    logs.sort(reverse=True)
    for v in logs:
        if w1 <= v <= w2:
            return True
    A = np.zeros(1); B = np.zeros(1)
    for i, v in enumerate(logs):
        if i % 2 == 0:
            A = np.concatenate((A, A + v)); A = A[A <= w2]
        else:
            B = np.concatenate((B, B + v)); B = B[B <= w2]
    B.sort()
    lo = np.searchsorted(B, w1 - A, side="left")
    hi = np.searchsorted(B, w2 - A, side="right")
    return bool(np.any(hi > lo))

def worker(args):
    x, keep_cut, n_target, seed, time_budget = args
    _init_sieve()
    rng = random.Random(seed)
    hits = tot = attempts = 0
    t0 = time.time()
    while tot < n_target and time.time() - t0 < time_budget:
        attempts += 1
        s = rng.randrange(1, x + 1)
        m = 1
        fac = []
        dead = False
        while s > 1:
            if is_prime_int(s):
                m *= s
                if m > x:
                    dead = True
                    break
                fac.append(s)
            s = rng.randrange(1, s + 1)
        if dead:
            continue
        if rng.randrange(x) >= m:
            continue
        if m <= keep_cut:
            continue
        tot += 1
        Lm = math.log(m)
        if subset_hit([math.log(p) for p in fac], Lm/2 - LN2, Lm/2 + LN2):
            hits += 1
    return hits, tot, attempts

def append_row(out_csv, tag, lnx, n, hits, att, secs):
    new = not os.path.exists(out_csv)
    with open(out_csv, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["tag","lnx","n","hits","rho","se","attempts","seconds"])
        rho = hits / n if n else float("nan")
        se = math.sqrt(rho * (1 - rho) / n) if n else float("nan")
        w.writerow([tag, lnx, n, hits, rho, se, att, round(secs, 1)])

def run_height(tag, N, min_hits, procs, hours, round_mins, seed_base,
               out_csv, hits_done, n_done, start_round=0):
    k = N // 2
    x = k * k
    lnx = math.log(x)
    keep_cut = int(x * math.exp(-1.0))
    t_start = time.time()
    rnd = start_round
    if n_done == 0:
        h0, t0n, a0 = worker((x, keep_cut, 10**9, seed_base, 30.0))
        if t0n:
            append_row(out_csv, tag, lnx, t0n, h0, a0, 30.0)
        hits_done, n_done = h0, t0n
        rate = max(t0n / 30.0, 1e-4)
        rho_g = max(hits_done / max(n_done, 1), 0.02)
        eta_h = max(min_hits - hits_done, 0) / rho_g / (rate * procs) / 3600
        print(f"{tag}: warm-up {t0n} samples ({rate:.2f}/s/core), "
              f"rho~{rho_g:.3f}, ETA {eta_h:.2f} h on {procs} cores",
              flush=True)
    else:
        print(f"{tag}: resuming with {hits_done} hits / {n_done} samples "
              f"banked", flush=True)
    while hits_done < min_hits and time.time() - t_start < hours * 3600:
        rnd += 1
        rho_g = max(hits_done / max(n_done, 1), 0.02)
        n_round = max(200, int((min_hits - hits_done) / rho_g))
        per = -(-n_round // procs)
        rbudget = min(round_mins * 60.0,
                      hours * 3600 - (time.time() - t_start))
        args = [(x, keep_cut, per, seed_base + 10**4 * rnd + i, rbudget)
                for i in range(procs)]
        r0 = time.time()
        with mp.Pool(procs) as pool:
            res = pool.map(worker, args)
        rh = sum(r[0] for r in res); rn = sum(r[1] for r in res)
        ra = sum(r[2] for r in res)
        if rn:
            append_row(out_csv, tag, lnx, rn, rh, ra, time.time() - r0)
        hits_done += rh; n_done += rn
        print(f"{tag} round {rnd}: +{rn} samples, banked "
              f"{hits_done}/{min_hits} hits "
              f"(rho={hits_done/max(n_done,1):.5f})", flush=True)
    rho = hits_done / n_done if n_done else float("nan")
    print(f"{tag}: done/paused at rho = {rho:.5f} "
          f"(n={n_done}, hits={hits_done})", flush=True)

def read_banked(out_csv):
    banked = {}
    if os.path.exists(out_csv):
        with open(out_csv) as f:
            for r in csv.DictReader(f):
                h, n, rc = banked.get(r["tag"], (0, 0, 0))
                banked[r["tag"]] = (h + int(r["hits"]), n + int(r["n"]),
                                    rc + 1)
    return banked

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-hits", type=int, default=400)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--hours", type=float, default=12.0)
    ap.add_argument("--round-mins", type=float, default=20.0)
    ap.add_argument("--out", default="mc_tier1.csv")
    ap.add_argument("--seed", type=int, default=20260822)
    ap.add_argument("--redo", action="store_true",
                    help="rerun heights even if banked hits >= min-hits")
    a = ap.parse_args()
    _init_sieve()
    banked = read_banked(a.out)
    heights = [("N=10^8", 10**8), ("2^64", 2**64), ("2^128", 2**128),
               ("2^256", 2**256), ("2^512", 2**512), ("2^1024", 2**1024)]
    for i, (tag, N) in enumerate(heights):
        h, n, rc = banked.get(tag, (0, 0, 0))
        target = h + a.min_hits if a.redo else a.min_hits
        if h >= target:
            print(f"{tag}: {h} hits banked >= {target}, skipping "
                  f"(--redo adds another --min-hits worth)", flush=True)
            continue
        run_height(tag, N, target, a.procs, a.hours, a.round_mins,
                   a.seed + 10**7 * i, a.out, h, n, start_round=rc)

if __name__ == "__main__":
    main()
