"""Bound and estimate of gamma - Delta from the volatility of the one-year real rate measured in advance
(Proposition 4: r^f_{t+1} = const - (gamma - Delta) * eps_t, eps_t the surprise in consumption growth of year t).

Data (FRED, in data/FRED/): GS1 (1-year Treasury constant maturity, monthly), MICH (University of Michigan median
expected inflation over the next 12 months, monthly, from 1978), A794RX0Q048SBEA (real PCE per capita, chained 2017
dollars, quarterly). Real rate measured in advance: GS1 - MICH, annual average of monthly values (and, as a check,
the January value, set right after the year whose consumption growth is used). Consumption growth: log change of the
annual average of quarterly real PCE per capita.

With i.i.d. growth (the model of Proposition 4):
  (1) sd(r^f) = (gamma - Delta) sigma, so gamma - Delta <= sd(r^f)/sigma (any other source of variation in the rate,
      e.g. monetary policy or slow-moving trends, only loosens this as an upper bound: it is an upper bound);
  (2) in first differences, sd(Delta r^f) = sqrt(2) (gamma - Delta) sigma, which removes slow trends;
  (3) the regression of r^f_{t+1} on g_t has slope -(gamma - Delta); in differences, Delta r^f_{t+1} on Delta g_t, too.

Writes out/data/real_rate_bound.csv (statistics) and out/data/real_rate_bound_data.csv (the annual series).

Superseded by bound_curve.py, which computes the same bound at two aspiration horizons, measures the
surprise against survey expectations rather than against the unconditional standard deviation of growth,
and produces Figure 2.  No number in the paper comes from this script; it is kept because it is the
simplest annual version of the calculation and uses FRED data only, with no survey workbooks, so it runs
on its own and serves as an independent check on the annual i.i.d.-sigma column of out/data/bound_curve.csv."""
import csv, os
import numpy as np
here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fred = os.path.join(here, "data", "FRED")
out = os.path.join(here, "out")
# Machine-readable output (CSV) goes in out/data; text reports stay in out/.
odata = os.path.join(out, "data"); os.makedirs(odata, exist_ok=True)

def load(name):
    """Read one FRED CSV and return {observation date as a string: value}, skipping missing values."""
    with open(os.path.join(fred, name + ".csv")) as f:
        r = csv.reader(f); next(r)
        return {d: float(v) for d, v in r if v not in ("", ".")}

gs1, mich, pce = load("GS1"), load("MICH"), load("A794RX0Q048SBEA")
# annual real rate: average of the twelve months, and the January value
rr = {d: gs1[d] - mich[d] for d in gs1 if d in mich}
years = sorted({d[:4] for d in rr})
r_avg = {y: np.mean([rr[d] for d in rr if d[:4] == y]) for y in years if sum(d[:4] == y for d in rr) == 12}
r_jan = {y: rr[y + "-01-01"] for y in years if y + "-01-01" in rr}
# annual consumption growth (percent)
cy = {}
for d, v in pce.items(): cy.setdefault(d[:4], []).append(v)
c_avg = {y: np.mean(v) for y, v in cy.items() if len(v) == 4}
g = {y: 100*np.log(c_avg[y]/c_avg[str(int(y) - 1)]) for y in c_avg if str(int(y) - 1) in c_avg}

def newey_west(X, e, lags):
    """Newey-West (Bartlett kernel) standard errors of the OLS coefficients of X on residuals e."""
    n = len(e); XtXi = np.linalg.inv(X.T @ X); u = X * e[:, None]; S = u.T @ u
    for l in range(1, lags + 1):
        w = 1 - l/(lags + 1); G = u[l:].T @ u[:-l]; S += w*(G + G.T)
    return np.sqrt(np.diag(XtXi @ S @ XtXi))

def ols(y, x, lags=2):
    """Least squares of y on a constant and x; returns the slope and its Newey-West standard error."""
    X = np.column_stack([np.ones(len(x)), x]); b = np.linalg.lstsq(X, y, rcond=None)[0]
    return b[1], newey_west(X, y - X @ b, lags)[1]

rows = []
for rname, r in (("annual average", r_avg), ("January", r_jan)):
    for excl in (False, True):
        T = [y for y in sorted(r) if str(int(y) - 1) in g and (not excl or y not in ("2020", "2021", "2022"))]
        rt = np.array([r[y] for y in T]); gt = np.array([g[str(int(y) - 1)] for y in T])   # r_{t+1} and g_t
        sig = np.std(gt, ddof=1); sdr = np.std(rt, ddof=1)
        drt, dgt = np.diff(rt), np.diff(gt)
        b_lev, se_lev = ols(rt, gt); b_dif, se_dif = ols(drt, dgt)
        rows.append(dict(rate=rname, sample=f"{T[0]}-{T[-1]}" + (" excl. 2020-22" if excl else ""), n=len(T),
                         sd_consumption_growth=sig, sd_real_rate=sdr, bound_levels=sdr/sig,
                         sd_change_real_rate=np.std(drt, ddof=1), bound_changes=np.std(drt, ddof=1)/(np.sqrt(2)*sig),
                         slope_levels=b_lev, se_levels=se_lev, implied_gap_levels=-b_lev,
                         slope_changes=b_dif, se_changes=se_dif, implied_gap_changes=-b_dif))
with open(os.path.join(odata, "real_rate_bound.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(os.path.join(odata, "real_rate_bound_data.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["year", "real_rate_avg", "real_rate_january", "consumption_growth"])
    for y in sorted(set(r_avg) | set(g)): w.writerow([y, r_avg.get(y, ""), r_jan.get(y, ""), g.get(y, "")])
for x in rows:
    print(f"{x['rate']:15s} {x['sample']:22s} n={x['n']:2d} sd(g)={x['sd_consumption_growth']:.2f} sd(r)={x['sd_real_rate']:.2f} "
          f"bound={x['bound_levels']:.2f} bound(changes)={x['bound_changes']:.2f} | slope levels {x['slope_levels']:+.2f} ({x['se_levels']:.2f}) "
          f"changes {x['slope_changes']:+.2f} ({x['se_changes']:.2f})")
