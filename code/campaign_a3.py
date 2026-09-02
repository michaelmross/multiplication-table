#!/usr/bin/env python3
"""
Campaign A3: exact-arithmetic rerun of the deep arm (2^10^17, 2^10^18).

Why. The longdouble window test blurs subset sums by ~0.05-0.5 at
L = 1.4e17..1.4e18, which a dose-response test showed shifts the measured
density by tens of percent -- voiding the original deep-arm cells.  This
engine removes all summation rounding from the window classification:

  * Each stick w (generated in longdouble exactly as before -- the generator
    was exonerated by the dl/l small-stick diagnostic) is split EXACTLY into
    (floor(w): integer, w - floor(w): float64 fraction; the split is exact
    in longdouble, and the one-time float64 cast of the fraction merely
    redefines the stick by <= 1e-17 -- a value perturbation, not an
    accumulating error).
  * Subset sums carry (int64 integer part, float64 fraction part) with
    carry-normalization; integer parts are EXACT (sums <= L < 2^62), and
    fraction parts (bounded by the stick count) are exact to ~1e-14.
  * Window comparisons are exact integer comparisons except within one unit
    of an edge, where they are float64-exact to ~1e-13.
  * Meet-in-the-middle range queries use a dense-rank composite key
    G = 2*rank(int part) + frac -- exact in float64 because ranks are small
    -- so the search is vectorized and rounding-free.
  Ceiling of this representation: L < 2^62 ~ 4.6e18.

MODES
  validate : mandatory gate.  Runs exact vs longdouble on IDENTICAL stick
             multisets: at L = 1.39e12 they must agree on >= 99.95% of
             samples (both engines are sound there); at L = 1.39e17 and
             1.39e18 the disagreement rate is REPORTED -- it is the direct
             measurement of the artifact that voided the original cells.
             Writes .a3_validated on success; 'run' refuses without it.
  run      : top-up 2^10^17 and 2^10^18 at c in {0.5, 0.75} to 20000 hits
             per cell.  Writes to mc_tier3.csv ONLY and counts banked hits
             ONLY from mc_tier3.csv -- the voided tier-2 rows for these
             tags are never counted or overwritten.  Kill/relaunch safe;
             seeds advance with banked rows.
  evaluate : frozen D1/D2/D3 (C_frozen = 0.0126, b_frozen = 0.359,
             unchanged) on the merged clean set: mc_tier2.csv MINUS tags
             {2^10^17, 2^10^18} UNION mc_tier3.csv.

Usage:
    python3 campaign_a3.py validate
    python3 campaign_a3.py run [--scale S] [--procs P] [--round-mins M]
    python3 campaign_a3.py evaluate
"""
import argparse, csv, math, os, random, time
import multiprocessing as mp
import numpy as np

LD = np.longdouble
LN2 = math.log(2.0)
DELTA = 1.0 - (1.0 + math.log(math.log(2.0))) / math.log(2.0)
C_FROZEN = 0.0126
B_FROZEN = 0.359
STAMP = ".a3_validated"

# ---------------- stick generation (unchanged, exonerated) ----------------
def gen_sticks(L, rng):
    rem = L
    s = []
    ln2 = LD(LN2)
    while rem >= ln2:
        w = rem * LD(rng.random())
        if w >= ln2:
            s.append(w)
        rem -= w
    return s

# ---------------- exact (int, frac) machinery ----------------
def split(sticks):
    ai = []; af = []
    for w in sticks:
        a = int(w)                     # floor for positive w
        ai.append(a); af.append(float(w - LD(a)))
    return ai, af

def norm_scalar(i, f):
    c = math.floor(f)
    return i + c, f - c

def val_le(i1, f1, i2, f2):
    """(i1+f1) <= (i2+f2), fracs in [0,1)."""
    if i1 != i2:
        return i1 < i2
    return f1 <= f2

def hit_exact(sticks, h):
    ai, af = split(sticks)
    Mi = sum(ai); Mf = math.fsum(af)
    W1i, W1f = norm_scalar(Mi // 2, (Mi % 2) * 0.5 + Mf / 2 - h)
    W2i, W2f = norm_scalar(Mi // 2, (Mi % 2) * 0.5 + Mf / 2 + h)
    keep = [(a, f) for a, f in zip(ai, af) if val_le(a, f, W2i, W2f)]
    if not keep:
        return False
    Ti, Tf = norm_scalar(sum(a for a, _ in keep),
                         math.fsum(f for _, f in keep))
    if not val_le(W1i, W1f, Ti, Tf):
        return False
    keep.sort(key=lambda p: p[0] + p[1], reverse=True)
    for a, f in keep:
        if val_le(W1i, W1f, a, f) and val_le(a, f, W2i, W2f):
            return True
    if len(keep) > 36:
        return dfs_exact(keep, (W1i, W1f), (W2i, W2f))
    return mim_exact(keep, (W1i, W1f), (W2i, W2f))

def dfs_exact(keep, W1, W2):
    K = len(keep)
    sufi = [0] * (K + 1); suff = [0.0] * (K + 1)
    for j in range(K - 1, -1, -1):
        sufi[j] = sufi[j + 1] + keep[j][0]
        suff[j] = suff[j + 1] + keep[j][1]
    stack = [(0, 0, 0.0)]
    while stack:
        j, si, sf = stack.pop()
        si, sf = norm_scalar(si, sf)
        if val_le(W1[0], W1[1], si, sf) and val_le(si, sf, W2[0], W2[1]):
            return True
        if j == K or not val_le(si, sf, W2[0], W2[1]):
            continue
        ti, tf = norm_scalar(si + sufi[j], sf + suff[j])
        if not val_le(W1[0], W1[1], ti, tf):
            continue
        stack.append((j + 1, si, sf))
        stack.append((j + 1, si + keep[j][0], sf + keep[j][1]))
    return False

def _build_half(items, W2i, W2f):
    Ai = np.zeros(1, dtype=np.int64); Af = np.zeros(1)
    for a, f in items:
        Ai = np.concatenate((Ai, Ai + a)); Af = np.concatenate((Af, Af + f))
        c = np.floor(Af); Ai = Ai + c.astype(np.int64); Af = Af - c
        m = (Ai < W2i) | ((Ai == W2i) & (Af <= W2f))
        Ai = Ai[m]; Af = Af[m]
    return Ai, Af

def mim_exact(keep, W1, W2):
    W1i, W1f = W1; W2i, W2f = W2
    Ai, Af = _build_half(keep[0::2], W2i, W2f)
    Bi, Bf = _build_half(keep[1::2], W2i, W2f)
    if len(Ai) == 0 or len(Bi) == 0:
        return False
    order = np.lexsort((Bf, Bi))
    Bi = Bi[order]; Bf = Bf[order]
    U = np.unique(Bi)
    ranks = np.searchsorted(U, Bi)
    G = 2.0 * ranks + Bf                       # exact: ranks small
    LOi = W1i - Ai; LOf = W1f - Af
    c = np.floor(LOf); LOi = LOi + c.astype(np.int64); LOf = LOf - c
    HIi = W2i - Ai; HIf = W2f - Af
    c = np.floor(HIf); HIi = HIi + c.astype(np.int64); HIf = HIf - c
    def count(Qi, Qf, right):
        r = np.searchsorted(U, Qi, side="left")
        inb = r < len(U)
        match = np.zeros(len(Qi), dtype=bool)
        match[inb] = U[r[inb]] == Qi[inb]
        if right:
            q = 2.0 * r + np.where(match, np.nextafter(Qf, 2.0), 0.0)
        else:
            q = 2.0 * r + np.where(match, Qf, 0.0)
        return np.searchsorted(G, q, side="left")
    lo = count(LOi, LOf, right=False)
    hi = count(HIi, HIf, right=True)
    return bool(np.any(hi > lo))

def one_exact(L, h, rng):
    return hit_exact(gen_sticks(L, rng), h)

def worker(args):
    L, h, n, seed, tcap = args
    rng = random.Random(seed)
    t0 = time.time(); hits = tot = 0
    while tot < n and time.time() - t0 < tcap:
        for _ in range(min(20, n - tot)):
            hits += one_exact(L, h, rng); tot += 1
    return hits, tot

# ---------------- validate ----------------
def pair_worker(args):
    """Paired exact-vs-longdouble comparison on identical stick streams.
    Returns (n, both_hit, exact_only, ld_only)."""
    from campaign_a import hit as hit_ld
    Lf, hf, n, seed = args
    rng = random.Random(seed)
    L = LD(Lf)
    both = e_only = l_only = 0
    for _ in range(n):
        s = gen_sticks(L, rng)
        M = LD(0.0)
        for v in s:
            M += v
        e = hit_exact(s, hf)
        l = hit_ld(list(s), M / 2 - LD(hf), M / 2 + LD(hf))
        both += (e and l); e_only += (e and not l); l_only += (l and not e)
    return n, both, e_only, l_only

def validate(procs):
    """Two jobs in one paired design.  (1) GATE at a clean height: the
    engines must agree on >= 99.95% of samples.  (2) ARBITER at the deep
    heights, run at 60k pairs per (L, c) so that flip rates of order 1e-4
    -- the size the blur hypothesis needs to explain the deep-arm anomaly
    -- are visible.  Flip DIRECTIONS are reported: ld_only > exact_only
    means longdouble over-counts hits (blur inflation), the reverse means
    suppression.  If flips are absent at the deep heights, longdouble is
    vindicated there and the deep-arm anomaly must be treated as possibly
    REAL, not instrumental."""
    ok = True
    jobs = [(1.39e12, math.log(2.0), 20000, 0.9995, "gate"),
            (1.39e17, math.log(2.0), 60000, None, "arbiter c=0.5"),
            (1.39e18, math.log(2.0), 60000, None, "arbiter c=0.5"),
            (1.39e17, math.log(4/3), 60000, None, "arbiter c=0.75"),
            (1.39e18, math.log(4/3), 60000, None, "arbiter c=0.75")]
    for Lf, hf, n, gate, label in jobs:
        per = -(-n // procs)
        args = [(Lf, hf, per, 20260830 + int(Lf % 997) * 1000
                 + int(hf * 100) + j) for j in range(procs)]
        t0 = time.time()
        with mp.Pool(procs) as pool:
            res = pool.map(pair_worker, args)
        nt = sum(r[0] for r in res); both = sum(r[1] for r in res)
        e1 = sum(r[2] for r in res); l1 = sum(r[3] for r in res)
        flips = e1 + l1
        rate = 1 - flips / nt
        print(f"[{label}] L={Lf:.2e}: agreement {rate*100:.4f}% "
              f"({flips} flips: exact-only {e1}, ld-only {l1})  "
              f"rho_exact={(both+e1)/nt:.5f} rho_ld={(both+l1)/nt:.5f}  "
              f"[{(time.time()-t0)/60:.1f} min]", flush=True)
        if gate is not None and rate < gate:
            ok = False
            print(f"  GATE FAILED (required {gate*100:.2f}%)")
    if ok:
        open(STAMP, "w").write(str(time.time()))
        print("validation stamp written; 'run' is unlocked")
        print("interpretation: flips >~ 10 with one-sided direction at a "
              "deep height = blur artifact confirmed there; flips ~ 0 = "
              "longdouble vindicated and the deep-arm anomaly is possibly "
              "real -- either way, proceed with 'run': the exact engine is "
              "correct in both worlds.")
    else:
        print("validation FAILED; fix before running")

# ---------------- run ----------------
def read_banked(path):
    banked = {}
    if os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                k = (r["tag"], float(r["c"]))
                h, n, rc = banked.get(k, (0.0, 0, 0))
                banked[k] = (h + float(r["rho"]) * int(r["n"]),
                             n + int(r["n"]), rc + 1)
    return banked

def run(a):
    if not os.path.exists(STAMP):
        raise SystemExit("run 'python3 campaign_a3.py validate' first")
    new = not os.path.exists("mc_tier3.csv")
    f = open("mc_tier3.csv", "a", newline=""); w = csv.writer(f)
    if new:
        w.writerow(["tag","c","L","t","rho","se","n","seconds"])
    for e in (17, 18):
        tag = f"2^10^{e}"
        lnN = 10.0**e * LN2
        for c in (0.5, 0.75):
            target = int(20000 * a.scale)
            Lf = 2 * lnN + 2 * math.log(c)
            L = LD(Lf); h = math.log(1.0 / c)
            hb, nb, rc = read_banked("mc_tier3.csv").get((tag, c),
                                                         (0.0, 0, 0))
            if hb >= target:
                print(f"{tag} c={c}: {hb:.0f} banked >= {target}, skip",
                      flush=True)
                continue
            print(f"{tag} c={c}: banked {hb:.0f}/{target}", flush=True)
            while hb < target:
                rc += 1
                rho_g = max(hb / max(nb, 1), 5e-4)
                per = -(-max(2000, int((target - hb) / rho_g)) // a.procs)
                args = [(L, h, per,
                         a.seed + hash((tag, c)) % 10**6 + 10**7 * rc + j,
                         a.round_mins * 60.0) for j in range(a.procs)]
                t0 = time.time()
                with mp.Pool(a.procs) as pool:
                    res = pool.map(worker, args)
                rh = sum(r[0] for r in res); rn = sum(r[1] for r in res)
                if rn == 0:
                    break
                rho = rh / rn
                w.writerow([tag, c, Lf, math.log(Lf), rho,
                            math.sqrt(rho * (1 - rho) / rn), rn,
                            round(time.time() - t0, 1)]); f.flush()
                hb += rh; nb += rn
                print(f"{tag} c={c} round: +{rn}, banked {hb:.0f}/{target} "
                      f"(rho={hb/nb:.5f})", flush=True)
    f.close()

# ---------------- evaluate ----------------
def evaluate():
    rows = [dict(r) for r in csv.DictReader(open("mc_tier2.csv"))
            if r["tag"] not in ("2^10^17", "2^10^18")]
    if os.path.exists("mc_tier3.csv"):
        rows += [dict(r) for r in csv.DictReader(open("mc_tier3.csv"))]
    acc = {}
    for r in rows:
        k = (r["tag"], float(r["c"]))
        n = int(r["n"]); h = float(r["rho"]) * n
        if k in acc:
            acc[k]["n"] += n; acc[k]["h"] += h
        else:
            acc[k] = {"t": float(r["t"]), "n": n, "h": h}
    pts = []
    for c in (0.5, 0.75):
        sub = sorted((v for (tag, cc), v in acc.items() if cc == c),
                     key=lambda v: v["t"])
        t = np.array([v["t"] for v in sub])
        rho = np.array([v["h"] / v["n"] for v in sub])
        sig = np.sqrt((1 - rho) / (rho * np.array([v["n"] for v in sub])))
        tm = (t[1:] + t[:-1]) / 2
        sl = np.diff(np.log(rho)) / np.diff(t)
        se = np.sqrt(sig[1:]**2 + sig[:-1]**2) / np.diff(t)
        for a_, b_, s_ in zip(tm, sl, se):
            pts.append((a_, DELTA - (-(b_ + 1.5 / a_)), s_))
    for lo, hi, tg in ((20, 28.5, .0015), (28, 35.5, .0011),
                       (36.5, 42.5, .0013)):
        q = [(a, d, s) for a, d, s in pts if lo <= a < hi]
        if not q:
            print(f"D1 [{lo},{hi}): no data"); continue
        t_, d_, s_ = map(np.array, zip(*q))
        wg = 1 / s_**2
        tb = float((t_ * wg).sum() / wg.sum())
        db = float((d_ * wg).sum() / wg.sum())
        sb = 1 / math.sqrt(wg.sum())
        print(f"D1 [{lo},{hi}): t_bar={tb:.1f} deficit={db:+.4f}+-{sb:.4f} "
              f"(target {tg}: {'MET' if sb <= tg else 'not yet'})")
        if lo == 36.5:
            print(f"D2: vs H_const ({C_FROZEN:+.4f}): "
                  f"z={(db-C_FROZEN)/sb:+.2f}; vs H_1/t "
                  f"({B_FROZEN/tb:+.4f}): z={(db-B_FROZEN/tb)/sb:+.2f}")
    q = [(a, d, s) for a, d, s in pts if a >= 20]
    t_, d_, s_ = map(np.array, zip(*q))
    X = np.vstack([np.ones_like(t_), 1 / t_]).T
    W = np.diag(1 / s_**2)
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ X.T @ W @ d_
    print(f"D3: C_inf = {beta[0]:+.4f} +- {math.sqrt(cov[0,0]):.4f}   "
          f"b = {beta[1]:+.3f} +- {math.sqrt(cov[1,1]):.3f}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["validate", "run", "evaluate"])
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--round-mins", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=20260831)
    a = ap.parse_args()
    if a.mode == "validate":
        validate(a.procs)
    elif a.mode == "run":
        run(a)
    else:
        evaluate()
