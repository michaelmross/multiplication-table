#!/usr/bin/env python3
"""
Strip-multiplicity census for J_n = [4n^2 - n, 4n^2 + n], N = 2n.

R_N(m) := #{ a <= N : a | m }   (multiplicity of m in the N-row strip of the
multiplication table, counting the trivial row a = 1).

Facts used:
  * Every composite m in J_n has a prime factor <= N, and no m in J_n is a
    product of two integers both > N (since (2n+1)^2 > 4n^2 + n).
    Hence R_N(m) = 1  <=>  m is prime.

PREREGISTERED PREDICTIONS (stated before the run):
  P1. Mean multiplicity:  mean_{m in J_n} R_N(m) = H_N + O(1/n) * N-scale slop,
      H_N = sum_{a<=N} 1/a.  Predicted |mean - H_N| = O(1) with small constant;
      we record the deviation.
  P2. PNT rate for the lower tail:  #{R_N = 1} / u' -> 1, where
      u' = (2n+1) / log(4n^2).
  P3. Independent-divisibility (Mertens) model overshoots by exactly 2*e^{-gamma}:
      M(n) := (2n+1) * prod_{p <= N} (1 - 1/p)  satisfies
      M(n) / #{R_N = 1}  ->  2*e^{-gamma} = 1.1229189...
      (the numerical face of the Mertens/PNT discrepancy; the model's failure
      at exactly this factor is the toy shadow of the parity obstruction).
  P4. Window typicality (H0): the distribution of R_N over the square-centered
      window J_n does not differ detectably from the same statistic over a
      generic shifted window S_n = [4n^2 + 10n, 4n^2 + 12n] of equal width at
      the same magnitude.  Comparison: total-variation distance of normalized
      histograms, and the R = 1 counts.  (Caveat recorded: R_N = 1 on S_n is
      "prime OR p^2 with p in (2n, 2n+3]", a set of size O(1); and one shifted
      window is a single sample, so this is a diagnostic, not a test with
      power.)

Anchors: n in {10^3, 10^4, 10^5, 10^6, 10^7}  (last one skipped if slow).
Outputs: per-anchor histogram CSVs, summary CSV, and stdout report.
"""
import numpy as np
import math
import time
import csv
import sys

GAMMA = 0.5772156649015328606
TWO_E_MINUS_GAMMA = 2.0 * math.exp(-GAMMA)

def strip_multiplicity(L, W, N):
    """R[i] = #{a <= N : a | (L+i)} for i in [0, W)."""
    R = np.ones(W, dtype=np.int32)
    for a in range(2, N + 1):
        start = (-L) % a
        R[start::a] += 1
    return R

def prime_sieve(N):
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]:
            s[p*p::p] = False
    return np.nonzero(s)[0]

def harmonic(N):
    # exact-enough harmonic number via summation in float64 (N <= 2e7 fine)
    return np.sum(1.0 / np.arange(1, N + 1, dtype=np.float64))

def tv_distance(h1, h2):
    """Total variation distance between two integer-count histograms."""
    kmax = max(len(h1), len(h2))
    a = np.zeros(kmax); b = np.zeros(kmax)
    a[:len(h1)] = h1; b[:len(h2)] = h2
    a = a / a.sum(); b = b / b.sum()
    return 0.5 * np.abs(a - b).sum()

def census(n):
    N = 2 * n
    W = 2 * n + 1
    L = 4 * n * n - n            # J_n
    Ls = 4 * n * n + 10 * n      # shifted generic window S_n, same width
    t0 = time.time()
    R = strip_multiplicity(L, W, N)
    Rs = strip_multiplicity(Ls, W, N)
    t1 = time.time()

    primes_leq_N = prime_sieve(N)
    log_mertens = np.sum(np.log1p(-1.0 / primes_leq_N.astype(np.float64)))
    mertens_density = math.exp(log_mertens)

    hist = np.bincount(R)
    hist_s = np.bincount(Rs)

    r1 = int(hist[1]) if len(hist) > 1 else 0
    r2 = int(hist[2]) if len(hist) > 2 else 0
    r1_s = int(hist_s[1]) if len(hist_s) > 1 else 0

    u_prime = W / math.log(4.0 * n * n)
    mertens_pred = W * mertens_density
    HN = harmonic(N)

    out = dict(
        n=n, N=N, W=W,
        primes=r1, r2=r2, primes_shifted=r1_s,
        u_prime=u_prime,
        pnt_ratio=r1 / u_prime,
        mertens_pred=mertens_pred,
        mertens_over_emp=mertens_pred / r1 if r1 else float('nan'),
        target_2eg=TWO_E_MINUS_GAMMA,
        mean_R=float(R.mean()), HN=HN, mean_minus_HN=float(R.mean() - HN),
        max_R=int(R.max()),
        tv_Jn_vs_shifted=tv_distance(hist, hist_s),
        seconds=t1 - t0,
    )
    return out, hist, hist_s

def main():
    anchors = [10**3, 10**4, 10**5, 10**6, 10**7]
    rows = []
    for n in anchors:
        # crude time guard: skip 10^7 if 10^6 took > 120 s
        if n == 10**7 and rows and rows[-1]['seconds'] > 120:
            print(f"skipping n={n}: previous anchor too slow", flush=True)
            continue
        out, hist, hist_s = census(n)
        rows.append(out)
        np.savetxt(f"/home/claude/hist_Jn_n{n}.csv",
                   np.column_stack([np.arange(len(hist)), hist]),
                   fmt="%d", delimiter=",", header="R,count", comments="")
        np.savetxt(f"/home/claude/hist_shifted_n{n}.csv",
                   np.column_stack([np.arange(len(hist_s)), hist_s]),
                   fmt="%d", delimiter=",", header="R,count", comments="")
        print(f"n={n:>9}  primes={out['primes']:>8}  "
              f"pnt_ratio={out['pnt_ratio']:.5f}  "
              f"mertens/emp={out['mertens_over_emp']:.5f} (target {TWO_E_MINUS_GAMMA:.5f})  "
              f"meanR-H_N={out['mean_minus_HN']:+.5f}  "
              f"TV(Jn,shift)={out['tv_Jn_vs_shifted']:.5f}  "
              f"[{out['seconds']:.1f}s]", flush=True)
    with open("/home/claude/summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

if __name__ == "__main__":
    main()
