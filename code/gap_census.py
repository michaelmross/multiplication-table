#!/usr/bin/env python3
"""
Distinct entries of the N x N multiplication table per square-gap.

For c in (0,1), k = floor(c*N), G_k = (k^2, (k+1)^2], define
    D(N, c) = #{ m in G_k : m = ab for some a, b <= N },
    rho(N, c) = D / |G_k|,  |G_k| = 2k + 1.
Membership test: m in the table  <=>  m has a divisor a with m/N <= a <= N
(then b = m/a <= N).  For m ~ c^2 N^2 the divisor window [m/N, N] has bounded
multiplicative ratio ~ 1/c^2, which is the regime of Ford's
H(x, y, Cy) ~ x / ((log y)^delta (loglog y)^{3/2}),  delta = 1 - (1+loglog2)/log2
                                                            = 0.086071...
Ford's theorem gives the AVERAGE of D over gaps at scale k ~ cN.  Whether
individual gaps track the average is unproven.  This census measures it.

PREREGISTERED PREDICTIONS (stated before the run):
  P1 (Ford-shape decay). With N0 = 10^3, the normalized density
      rho(N, c) / rho(N0, c)  tracks
      F(N) = (log N / log N0)^{-delta} * (loglog N / loglog N0)^{-3/2}
      to within 20% relative at every anchor and every c.
      Numerically F(10^8) = 0.497: the density should roughly HALVE.
      Discriminants: constant-density null predicts 1.00 (fails by ~2x);
      pure (log N)^{-delta} predicts 0.92 (fails by ~1.8x at the far anchor).
      NOTE (honesty): the decay at these heights is dominated by the
      (loglog)^{3/2} factor; delta itself contributes only ~8% over five
      decades.  This census can confirm or refute the Ford SHAPE; it cannot
      measure delta to better than roughly +/- 0.1 (see P3).
  P2 (Universality in c). The normalized decay curves for c = 0.5, 0.75, 0.9
      coincide within their mutual 20% bands (the exponents are independent
      of the window ratio 1/c^2; only Ford's constant depends on c).
      Level prediction: rho increases as c decreases (wider divisor window).
  P3 (Slope diagnostic, weak). Fitting
      log rho = A - d * log log N - (3/2) log(loglog N)   (3/2 held fixed)
      yields d in (0, 0.25) for every c.  Recorded as a diagnostic, not a
      measurement of delta, per the note in P1.
  P4 (Per-gap concentration). Over the ensemble of 10 consecutive gaps
      k, k+1, ..., k+9 (N <= 10^7; 3 gaps at N = 10^8 under the time guard),
      the coefficient of variation of D is < 10% at every anchor and is
      nonincreasing in N up to noise -- the empirical face of the unproven
      "almost all gaps" statement.

Anchors: N = 10^3, 10^4, 10^5, 10^6, 10^7, 10^8 (skip 10^8 if 10^7 slow).
Correctness: brute-force table construction cross-check at N in {97, 100, 211,
300} for all three c before the main run; abort on mismatch.
Exact integer arithmetic throughout the sieve (int64 with overflow headroom:
max products ~ 10^16 << 2^63).
"""
import numpy as np
import math
import time
import csv

DELTA = 1.0 - (1.0 + math.log(math.log(2.0))) / math.log(2.0)   # 0.086071...
CS = [0.5, 0.75, 0.9]
N0 = 10**3

def distinct_in_gap(N, k, chunk=5_000_000):
    """#{m in (k^2, (k+1)^2] : m = ab, a,b <= N}, exact, vectorized."""
    lo = k * k + 1                 # inclusive
    hi = (k + 1) * (k + 1)         # inclusive
    glen = hi - lo + 1
    hit = np.zeros(glen, dtype=bool)
    a_min = max(1, -(-lo // N))    # ceil(lo/N)
    a_lo = a_min
    while a_lo <= N:
        a_hi = min(N, a_lo + chunk - 1)
        a = np.arange(a_lo, a_hi + 1, dtype=np.int64)
        m = ((lo + a - 1) // a) * a          # first multiple >= lo
        cap = np.minimum(np.int64(hi), a * np.int64(N))
        while True:
            ok = m <= cap
            if not ok.any():
                break
            hit[m[ok] - lo] = True
            m = m + a
        a_lo = a_hi + 1
    return int(hit.sum()), glen

def brute_check():
    for N in (97, 100, 211, 300):
        prods = set()
        for a in range(1, N + 1):
            for b in range(a, N + 1):
                prods.add(a * b)
        for c in CS:
            k = int(c * N)
            lo, hi = k * k + 1, (k + 1) * (k + 1)
            truth = sum(1 for m in range(lo, hi + 1) if m in prods)
            got, _ = distinct_in_gap(N, k)
            assert got == truth, (N, c, got, truth)
    print("brute-force cross-check passed (N=97,100,211,300; all c)", flush=True)

def ford_factor(N):
    return ((math.log(N) / math.log(N0)) ** (-DELTA) *
            (math.log(math.log(N)) / math.log(math.log(N0))) ** (-1.5))

def main():
    brute_check()
    anchors = [10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
    rows = []
    slow = False
    for N in anchors:
        if N == 10**8 and slow:
            print("skipping 10^8: time guard", flush=True)
            continue
        for c in CS:
            k0 = int(c * N)
            n_gaps = 10 if N <= 10**7 else 3
            t0 = time.time()
            Ds, glens = [], []
            for j in range(n_gaps):
                D, glen = distinct_in_gap(N, k0 + j)
                Ds.append(D); glens.append(glen)
            dt = time.time() - t0
            Ds = np.array(Ds, dtype=float)
            rhos = Ds / np.array(glens, dtype=float)
            rho = rhos.mean()
            cv = Ds.std(ddof=1) / Ds.mean() if n_gaps > 1 else float('nan')
            rows.append(dict(N=N, c=c, k=k0, n_gaps=n_gaps,
                             D_gap0=int(Ds[0]), rho=rho, cv=cv,
                             ford_factor=ford_factor(N), seconds=dt))
            print(f"N=10^{int(math.log10(N))} c={c:.2f}  rho={rho:.5f}  "
                  f"CV={cv*100:.2f}%  [{dt:.1f}s]", flush=True)
        if rows[-1]['seconds'] * (10**8 // N) > 600 and N < 10**8:
            slow = (N == 10**7 and rows[-1]['seconds'] > 60)
    # normalized decay + slope diagnostic
    print("\nP1/P2: normalized decay rho(N)/rho(10^3) vs Ford factor")
    base = {c: next(r['rho'] for r in rows if r['N'] == N0 and r['c'] == c)
            for c in CS}
    for r in rows:
        r['norm'] = r['rho'] / base[r['c']]
        r['norm_over_ford'] = r['norm'] / r['ford_factor']
    for c in CS:
        line = "  ".join(f"{r['norm']:.4f}/{r['ford_factor']:.4f}"
                         for r in rows if r['c'] == c)
        print(f"c={c:.2f}:  {line}")
    print("\nP3: fitted d (3/2 fixed), model log rho = A - d loglogN - 1.5 log(loglogN)")
    for c in CS:
        sub = [r for r in rows if r['c'] == c]
        x = np.array([math.log(math.log(r['N'])) for r in sub])
        y = np.array([math.log(r['rho']) + 1.5 * math.log(math.log(math.log(r['N'])))
                      for r in sub])
        d = -np.polyfit(x, y, 1)[0]
        print(f"c={c:.2f}: d_hat = {d:+.4f}   (delta = {DELTA:.4f})")
    with open("/home/claude/gap_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

if __name__ == "__main__":
    main()
