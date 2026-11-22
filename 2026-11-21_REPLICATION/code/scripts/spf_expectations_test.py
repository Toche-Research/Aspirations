"""Euler equation in the forward direction, with survey expectations of consumption growth.

From Proposition (Time-varying moments), solving the risk-free-rate equation for expected growth:
    aspirations : mu_t = (r_t + log beta)/Delta + k (g_t - mu_{t-1}),   k = (gamma-Delta)/Delta
    catching up : mu_t = (r_t + log beta)/gamma + chi g_t,              chi = (gamma-Delta)/gamma
    power       : mu_t = (r_t + log beta)/gamma
(variance terms omitted). Regression  mu_t = a + b_r r_t + b_g g_t + b_m mu_{t-1} + e_t:
    aspirations b_m = -b_g < 0;  catching up b_m = 0 < b_g;  power b_g = b_m = 0;
    sticky (partially updated) expectations, a non-preference alternative, b_m > 0.
Two quarters ahead, E_t g_{t+2} = (E_t r_{t+1} + log beta)/Delta under aspirations: no loading on (g_t, mu_{t-1});
under catching up it loads chi^2 on g_t.
Same data and timing as spf_aspirations_test.py (see its header for the required data layout).
Newey-West standard errors, 4 lags.

Feeds the 'Euler-equation tests' paragraph of Section 7: the forward regression of expected growth on the
real rate, the latest growth and the previous expectation, whose positive weight on the previous
expectation tells against both aspirations and catching up.  Writes out/data/spf_expectations_results.csv.
"""
import os
import numpy as np
import pandas as pd
from aspirations import ODATA
from aspirations.data import load
from aspirations.regress import ols


def run(d, label, y):
    """Regress the expectation y on the real rate, the latest growth and the previous expectation.

    Returns one row of results, including sum_g_mulag with its p-value (aspirations set the sum to zero,
    with a negative weight on the previous expectation) and p_both0, the joint test that both weights are
    zero, which power utility implies.
    """
    m, s = ols(d, y, ["r", "g", "mu_lag"])
    est, p_sum = m.lin([{"g": 1, "mu_lag": 1}])
    _, p_zero = m.lin([{"g": 1}, {"mu_lag": 1}])
    return {
        "spec": label, "dep": y, "N": m.nobs, "start": str(s.index.min()), "end": str(s.index.max()),
        "b_r": m.params["r"], "se_r": m.bse["r"],
        "b_g": m.params["g"], "se_g": m.bse["g"],
        "b_mulag": m.params["mu_lag"], "se_mulag": m.bse["mu_lag"],
        "sum_g_mulag": est[0], "p_sum0": p_sum, "p_both0": p_zero, "R2": m.rsquared,
    }


def main():
    """Run the six specifications, write out/data/spf_expectations_results.csv and print the table."""
    d = load()
    covid = (d.index >= pd.Period("2020Q1", freq="Q")) & (d.index <= pd.Period("2021Q2", freq="Q"))
    nc = d.copy()
    nc.loc[covid, ["mu", "mu_next", "g"]] = np.nan
    nc["mu_lag"] = nc["mu"].shift(1)
    pre = d.loc[:pd.Period("2007Q4", freq="Q")]
    rows = []
    for y, lab in [("mu", "one quarter ahead"), ("mu_next", "two quarters ahead")]:
        rows += [run(d, f"{lab}, full", y), run(nc, f"{lab}, excl. 2020Q1-2021Q2", y), run(pre, f"{lab}, to 2007Q4", y)]
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ODATA, "spf_expectations_results.csv"), index=False)
    pd.set_option("display.width", 250)
    print(res.round(3).to_string(index=False))
    return res


if __name__ == "__main__":
    main()
