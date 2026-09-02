#!/usr/bin/env python3
"""
Analysis for the per-gap density campaign.  Reads (whichever exist):
    gap_summary.csv   -- exact census (N <= 10^8)
    mc_tier1.csv      -- exact-integer Kalai MC (c = 1/2)
    mc_tier2.csv      -- continuum GEM model MC (c in {0.5, 0.75, 0.9})
Emits verdicts for the preregistered predictions and three charts.

Evaluation criteria (fixed in advance; see README_mc.md):
  P1  Tier-1 at N=10^8 vs exact rho = 0.25276 within max(2 SE, 1%).
  P2  Tier-2 vs exact at overlap heights: |model/exact - 1| <= 20%.
  Seam  Tier-2 vs Tier-1 at 2^64, 2^256, 2^1024 (c=0.5): offset reported.
  P3  Tier-2 local slopes vs Ford -delta - 1.5/t within 0.05 for t_mid >= 8;
      also report first t from which agreement holds and persists.
  P4  pairwise c-agreement of local slopes within 0.03 for t_mid >= 8.
  P5  fit -d - 1.5/t on t_mid in [10, 28]: d = 0.086 +/- 0.03.
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DELTA = 1.0 - (1.0 + math.log(math.log(2.0))) / math.log(2.0)
EXACT = {(10**4, 0.5): 0.31163, (10**4, 0.75): 0.20239, (10**4, 0.9): 0.09285,
         (10**6, 0.5): 0.27495, (10**6, 0.75): 0.18258, (10**6, 0.9): 0.08700,
         (10**8, 0.5): 0.25276, (10**8, 0.75): 0.17017, (10**8, 0.9): 0.08280}
TAG2N = {"N=10^4": 10**4, "N=10^6": 10**6, "N=10^8": 10**8}

def read(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [dict(r) for r in csv.DictReader(f)]

def pool_rows(rows, keyf):
    """Merge duplicate rows (resume rounds) by sample count."""
    acc = {}
    for r in rows:
        k = keyf(r)
        n = int(r["n"]); hits = r["rho"] * n
        if k in acc:
            acc[k]["n"] += n; acc[k]["hits"] += hits
        else:
            acc[k] = {**r, "n": n, "hits": hits}
    out = []
    for r in acc.values():
        if r["n"] == 0:
            continue
        r["rho"] = r["hits"] / r["n"]
        r["se"] = math.sqrt(max(r["rho"] * (1 - r["rho"]), 1e-12) / r["n"])
        out.append(r)
    return out

def main():
    t2 = read("mc_tier2.csv")
    for r in t2:
        r["c"] = float(r["c"]); r["t"] = float(r["t"])
        r["rho"] = float(r["rho"]); r["se"] = float(r["se"])
    t2 = pool_rows(t2, lambda r: (r["tag"], r["c"]))
    t1 = read("mc_tier1.csv")
    for r in t1:
        r["rho"] = float(r["rho"]); r["se"] = float(r["se"])
        r["lnx"] = float(r["lnx"])
    t1 = pool_rows(t1, lambda r: r["tag"])

    print("=== P1 (Tier-1 harness) ===")
    v = [r for r in t1 if r["tag"] == "N=10^8"]
    if v:
        r = v[0]; dev = abs(r["rho"] - 0.25276)
        tol = max(2*r["se"], 0.01*0.25276)
        print(f"  rho={r['rho']:.5f} vs exact 0.25276, |dev|={dev:.5f}, "
              f"tol={tol:.5f}: {'PASS' if dev <= tol else 'FAIL'}")
    else:
        print("  not run")

    print("=== P2 (model vs exact at overlap) ===")
    worst = 0.0
    for r in t2:
        if r["tag"] in TAG2N:
            ex = EXACT[(TAG2N[r["tag"]], r["c"])]
            off = r["rho"]/ex - 1
            worst = max(worst, abs(off))
            print(f"  {r['tag']} c={r['c']}: model/exact-1 = {off:+.3f}")
    if t2:
        print(f"  worst |offset| = {worst:.3f}: "
              f"{'PASS' if worst <= 0.20 else 'FAIL'} (bound 20%)")

    print("=== Seam (Tier-2 vs Tier-1, c=0.5) ===")
    for tag in ("2^64", "2^128", "2^256", "2^512", "2^1024"):
        a = [r for r in t1 if r["tag"] == tag]
        b = [r for r in t2 if r["tag"] == tag and r["c"] == 0.5]
        if a and b:
            off = b[0]["rho"]/a[0]["rho"] - 1
            print(f"  {tag}: model/integer-1 = {off:+.3f} "
                  f"(+- {b[0]['se']/a[0]['rho'] + a[0]['se']*b[0]['rho']/a[0]['rho']**2:.3f})")

    print("=== P3/P4/P5 (Tier-2 slopes) ===")
    slopes = {}
    for c in (0.5, 0.75, 0.9):
        sub = sorted((r for r in t2 if r["c"] == c), key=lambda r: r["t"])
        # average duplicate heights (re-runs) by weight
        t = np.array([r["t"] for r in sub])
        lr = np.log(np.array([r["rho"] for r in sub]))
        sig = np.array([r["se"]/r["rho"] for r in sub])
        tm = (t[1:]+t[:-1])/2
        sl = np.diff(lr)/np.diff(t)
        sse = np.sqrt(sig[1:]**2 + sig[:-1]**2)/np.diff(t)
        slopes[c] = (tm, sl, sse)
    p3_ok, first_t = True, None
    for c in (0.5, 0.75, 0.9):
        tm, sl, sse = slopes[c]
        for a, b, s in zip(tm, sl, sse):
            ford = -DELTA - 1.5/a
            d = b - ford
            mark = "OK " if abs(d) <= 0.05 else "DEV"
            if a >= 8 and abs(d) > 0.05:
                p3_ok = False
            print(f"  c={c} t={a:6.2f} slope={b:+.4f}(+-{s:.3f}) "
                  f"ford={ford:+.4f} dev={d:+.4f} {mark}")
    print(f"  P3 (all t>=8 within 0.05): {'PASS' if p3_ok else 'FAIL'}")
    tm0 = slopes[0.5][0]
    p4_ok = True
    for i, a in enumerate(tm0):
        if a < 8: continue
        vals = [slopes[c][1][i] for c in (0.5, 0.75, 0.9)
                if i < len(slopes[c][1])]
        if max(vals) - min(vals) > 0.03:
            p4_ok = False
    print(f"  P4 (pairwise c-agreement 0.03, t>=8): "
          f"{'PASS' if p4_ok else 'FAIL'}")
    pooled = []
    for c in (0.5, 0.75, 0.9):
        tm, sl, sse = slopes[c]
        m = (tm >= 10) & (tm <= 28)
        if m.any():
            d = -np.average(sl[m] + 1.5/tm[m], weights=1/sse[m]**2)
            pooled.extend(list(zip(sl[m], tm[m], sse[m])))
            print(f"  P5 c={c}: d_hat = {d:+.4f}")
    if pooled:
        sl_, tm_, se_ = map(np.array, zip(*pooled))
        d = -np.average(sl_ + 1.5/tm_, weights=1/se_**2)
        sed = 1/math.sqrt(np.sum(1/se_**2))
        ok = abs(d - DELTA) <= 0.03
        print(f"  P5 pooled: d_hat = {d:+.4f} +- {sed:.4f} vs delta {DELTA:.4f}: "
              f"{'PASS' if ok else 'FAIL'}")

    # charts
    cols = {0.5: "#1f77b4", 0.75: "#d62728", 0.9: "#2ca02c"}
    fig, ax = plt.subplots(figsize=(7.5, 5))
    tg = np.linspace(2.5, 29, 200)
    ax.plot(tg, -DELTA - 1.5/tg, "k--", label=r"Ford $-\delta-\frac{3}{2}/t$")
    ax.axhline(-DELTA, color="crimson", lw=.8, ls=":",
               label=r"$-\delta=-0.0861$")
    exact_t = [2.40, 2.86, 3.27, 3.63]
    exact_sl = {0.5: [-0.36, -0.32, -0.30, -0.29],
                0.75: [-0.28, -0.27, -0.25, -0.24],
                0.9: [-0.13, -0.15, -0.16, -0.17]}
    for c in (0.5, 0.75, 0.9):
        tm, sl, sse = slopes[c]
        ax.errorbar(tm, sl, yerr=sse, fmt="o-", ms=4, color=cols[c],
                    label=f"model, c={c}", capsize=2, lw=1)
        ax.plot(exact_t, exact_sl[c], "s", ms=5, mfc="none", color=cols[c])
    ax.plot([], [], "s", mfc="none", color="gray", label="exact census")
    ax.set_xlabel(r"$t = \log\log x$"); ax.set_ylabel(r"$d\log\rho/dt$")
    ax.set_title("Local decay rate of per-gap table density, to $N=2^{10^{11}}$")
    ax.legend(fontsize=9); fig.tight_layout()
    fig.savefig("mc_chart_slopes.png", dpi=150)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for c in (0.5, 0.75, 0.9):
        sub = sorted((r for r in t2 if r["c"] == c), key=lambda r: r["t"])
        ax.semilogy([r["t"] for r in sub], [r["rho"] for r in sub], "o-",
                    ms=4, color=cols[c], label=f"model, c={c}")
    for (N, c), v in EXACT.items():
        ax.semilogy([math.log(2*math.log(c*N))], [v], "s", ms=5,
                    mfc="none", color=cols[c])
    for r in t1:
        ax.semilogy([math.log(r["lnx"])], [r["rho"]], "D", ms=5,
                    color="black", mfc="none")
    ax.set_xlabel(r"$t=\log\log x$"); ax.set_ylabel(r"$\rho$")
    ax.set_title("Per-gap density across 23 e-folds of $\\log\\log$")
    ax.legend(fontsize=9); fig.tight_layout()
    fig.savefig("mc_chart_density.png", dpi=150)
    print("charts written: mc_chart_slopes.png, mc_chart_density.png")

if __name__ == "__main__":
    main()
