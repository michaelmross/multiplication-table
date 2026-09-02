#!/usr/bin/env python3
"""
Ballot skeleton for Ford's exponent: exact companion to the MC campaign.

In loglog coordinates the prime factors of a uniform integer form a
unit-rate Poisson process on [0, t], and Ford's mechanism reduces the
divisor-in-window density to the process staying above the line of slope
1/ln2 (subset-sum count 2^j keeping up with span).  The survival
probability S(t) is a killed negative-drift walk: slack -> slack +
Pois(ln2) - 1, killed below 0, started from Pois(theta) - 1.  Exact DP,
positive arithmetic, O(t^2).  S(t) ~ C e^{-delta t} t^{-3/2} with
delta = 1 - (1+ln ln 2)/ln 2 (large-deviation rate at a = 1/ln2) and the
universal -3/2 of conditioned negative-drift walks.

Validated against Monte Carlo (z = -0.33, +0.71 at t ~ 8, 12) and slack-cap
insensitive to machine precision.  Results (theta in [0.5, 2] as the
boundary-offset systematic):
  deficit(t~16) = 0.025-0.029        [full model measured: ~0.040]
  deficit(t~31) = 0.0091-0.0109      [measured: +0.0110 +- 0.0012]
  deficit(t~40) = 0.0063-0.0072      [measured: +0.0090 +- 0.0022]
  effective 1/t coefficient on t in [28,42]: c = 0.26-0.31
                                     [measured: b = 0.333 +- 0.032]
  and c slides to 0.21-0.24 on t in [45,62]: PREDICTION for any future
  deeper campaign -- the measured b should decline to ~0.22 by t ~ 50.
Usage: python3 ballot_skeleton.py [theta] [J]
"""
import math, sys
import numpy as np
LN2 = math.log(2.0)
DELTA = 1.0 - (1.0 + math.log(LN2)) / LN2

def survival_logs(theta=1.0, J=90, CAP=800):
    def pois(lam, kmax):
        p = np.zeros(kmax + 1); p[0] = math.exp(-lam)
        for k in range(1, kmax + 1):
            p[k] = p[k-1] * lam / k
        return p
    inc = pois(LN2, 200); init = pois(theta, CAP + 1)
    v = np.concatenate((init[1:CAP+1], [0.0]))
    lnS = math.log(v.sum()); out = [lnS]
    for _ in range(J):
        v = v / v.sum()
        w = np.convolve(v, inc)[:CAP+200][1:]
        lnS += math.log(w[:CAP].sum()); out.append(lnS)
        v = np.zeros(CAP + 1); v[:CAP] = w[:CAP]
    return np.array(out)

if __name__ == "__main__":
    theta = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    J = int(sys.argv[2]) if len(sys.argv) > 2 else 90
    lnS = survival_logs(theta, J)
    t = theta + LN2 * np.arange(J + 1)
    tm = (t[1:] + t[:-1]) / 2
    deficit = DELTA - (-(np.diff(lnS) / LN2) - 1.5 / tm)
    for a, d in zip(tm, deficit):
        print(f"t={a:6.2f}  deficit={d:+.5f}  c_eff={d*a:+.4f}")
