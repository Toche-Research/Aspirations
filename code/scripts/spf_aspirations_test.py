"""Test of Proposition (Time-varying moments) with survey expectations.

Model (isoelastic example, time-varying moments):
    r^f_{t+1} = const + Delta*mu_t + (gamma-Delta)*mu_{t-1} - (gamma-Delta)*g_t + variance terms
Regression:  r_t = a + b1*mu_t + b2*mu_{t-1} + b3*g_t + e_t
    aspirations : b2 = gamma-Delta > 0 and b2 + b3 = 0
    catching up : b2 = 0, b3 = -(gamma-Delta)
    power       : b2 = b3 = 0

Data (all annualized percentage points), Survey of Professional Forecasters (Philadelphia Fed), FRED:
    mu_t     = mean SPF forecast of real consumption growth in the survey quarter s (DRCONSUM2(s))
    mu_{t-1} = the same forecast in the previous survey (DRCONSUM2(s-1))
    g_t      = first-release real consumption growth of quarter s-1 (known to panelists at survey s)
    r_t      = 3-month T-bill (TB3MS, middle month of quarter s, the survey month)
               minus expected GDP-deflator inflation over the bill's life (mean of DPGDP2(s), DPGDP3(s))
Standard errors: Newey-West (HAC, Bartlett kernel), 4 lags; normal and chi-square p-values.  Levels and first differences.

Feeds the 'Euler-equation tests' paragraph of Section 7: the eleven specifications reported there, and the
coefficients and standard errors quoted from them.  Writes out/data/spf_test_results.csv (one row per
specification) and out/data/spf_test_data.csv (the assembled quarterly series).

Data layout.  The three Survey of Professional Forecasters workbooks are read from code/data/RCONSUM/ and the
deflator workbook from code/data/RCONSUM/PGDP/.  Both directories ship with the repository, so the script runs
from a clean checkout with no preparation.  bound_curve.py, spf_expectations_test.py and spf_pinned_rate.py
import the data root and the regression machinery from here.

This module also holds the regression machinery (Newey-West covariance, the Fit class and its Wald test)
used by spf_expectations_test.py, spf_pinned_rate.py and bound_curve.py: statsmodels is not a dependency.
"""
import os
import warnings

import numpy as np
import pandas as pd

from aspirations import ODATA
from aspirations.data import load, q_index
from aspirations.regress import Fit, ols

warnings.filterwarnings("ignore")

def run(d, label, xs=("mu", "mu_lag", "g"), extra=(), diff=False):
    """Run one specification of the real-rate regression and return its row of results.

    The first three regressors are expected growth, the previous expectation and the latest growth, whose
    coefficients the model restricts.  `extra` adds controls, `diff` takes first differences of everything.
    Reported alongside the coefficients: the sum of the second and third coefficients, which aspirations
    set to zero (p_sum0), the p-value of the previous expectation alone (p_mulag0), and the joint test that
    both are zero, which power utility implies (p_power).
    """
    xs = list(xs) + list(extra)
    dd = d.copy()
    y = "r"
    if diff:
        dd = dd[["r"] + xs].diff()
    m, s = ols(dd, y, xs)
    b_mu, b_lag, b_g = (m.params[x] for x in xs[:3])
    est_sum, p_sum = m.lin([{xs[1]: 1, xs[2]: 1}])
    _, p_pow = m.lin([{xs[1]: 1}, {xs[2]: 1}])
    row = {
        "spec": label, "N": int(m.nobs),
        "start": str(s.index.min()), "end": str(s.index.max()),
        "b_mu": b_mu, "se_mu": m.bse[xs[0]],
        "b_mulag": b_lag, "se_mulag": m.bse[xs[1]],
        "b_g": b_g, "se_g": m.bse[xs[2]],
        "sum_lag_g": est_sum[0], "p_sum0": p_sum,
        "p_mulag0": m.pvalues[xs[1]], "p_power": p_pow, "R2": m.rsquared,
    }
    for x in xs[3:]:
        row[f"b_{x}"] = m.params[x]
        row[f"se_{x}"] = m.bse[x]
    return row


def main():
    """Run the eleven specifications of Section 7, write the two CSV files and print the table."""
    d = load()
    d.to_csv(os.path.join(ODATA, "spf_test_data.csv"))
    pre = d.loc[:pd.Period("2007Q4", freq="Q")]
    nozlb = d[(d.tbill > 0.25)]
    covid = (d.index >= pd.Period("2020Q1", freq="Q")) & (d.index <= pd.Period("2021Q2", freq="Q"))
    nocovid = d.copy()
    nocovid.loc[covid, ["mu", "g"]] = np.nan          # drops 2020Q1-2021Q2 and the lags they feed
    rows = [
        run(d, "levels, full"),
        run(d, "differences, full", diff=True),
        run(pre, "levels, to 2007Q4"),
        run(pre, "differences, to 2007Q4", diff=True),
        run(nozlb, "levels, T-bill > 0.25"),
        run(nocovid.assign(mu_lag=nocovid.mu.shift(1)), "levels, excl. 2020Q1-2021Q2"),
        run(nocovid.assign(mu_lag=nocovid.mu.shift(1)), "differences, excl. 2020Q1-2021Q2", diff=True),
        run(d, "levels, mu = next quarter", xs=("mu_next", "mu_next_lag", "g")),
        run(d, "differences, mu = next quarter", xs=("mu_next", "mu_next_lag", "g"), diff=True),
        run(d, "levels, dispersion controls", extra=("disp", "disp_lag")),
        run(d, "differences, dispersion controls", extra=("disp", "disp_lag"), diff=True),
    ]
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ODATA, "spf_test_results.csv"), index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(res.round(3).to_string(index=False))
    print("g filled from second release in:", list(d.index[d.g_filled.fillna(False).astype(bool)].astype(str)))
    print("sample with all variables:", d[["r", "mu", "mu_lag", "g"]].dropna().index[[0, -1]])
    return res


if __name__ == "__main__":
    main()
