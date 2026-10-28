"""The Euler equation of the aspiration model read the other way: expected consumption growth given the real rate.

Proposition (time-varying moments), part (ii), with constant risk: the risk-free rate satisfies
    r = const + Delta*mu_t - (gamma - Delta)*(g_t - mu_{t-1}),
for ANY conditionally normal consumption process the consumer chooses (it is a restriction from the Euler equation,
as in Hall 1978). Solving for expected growth:
    mu_t = const + (1/Delta)*r + theta*eps_t,   theta = (gamma - Delta)/Delta,   eps_t = g_t - mu_{t-1}.
Power utility (Delta = gamma): theta = 0, Hall's random walk with drift. Aspirations: theta > 0, expected growth
rises after a positive surprise (momentum). Writing eps_t = g_t - mu_{t-1}, the coefficient on g_t is theta and the
coefficient on mu_{t-1} is -theta; sticky expectations (Carroll et al. 2020) instead put a positive weight on mu_{t-1}.
With time-varying risk, a term in sigma_t^2 is added (proxy: VIX^2).

Data (Survey of Professional Forecasters, Philadelphia Fed; FRED), quarterly, as in spf_aspirations_test.py:
    mu_t      = DRCONSUM2(s), the mean forecast of real consumption growth in the survey quarter s
    mu_{t-1}  = DRCONSUM2(s-1)
    g_t       = first-release growth of quarter s-1 (known at survey s)
    r         = 3-month T-bill (middle month of quarter s) minus expected GDP-deflator inflation
    sigma^2   = average of (VIX/100)^2 over quarter s (VXO before 1990)
Samples: all; excluding 2020-21; 1981-2007; excluding the zero-lower-bound quarters (2009Q1-2015Q4, 2020Q2-2022Q1).
Standard errors: Newey-West, 4 lags.

Feeds the 'Euler-equation tests' paragraph of Section 7: the coefficients on the real rate, on the surprise
entered directly, and on the surprise with the squared VIX added as a measure of risk.  Writes
out/data/spf_pinned_rate_results.csv.  See the header of spf_aspirations_test.py for the required data layout."""
import os
import numpy as np
import pandas as pd
from aspirations import DATA, ODATA
from aspirations.data import load
from aspirations.regress import Fit


def vix_quarterly():
    """Quarterly mean of the squared VIX, in percentage points squared divided by 100.

    The VXO is spliced on before the VIX begins, so the series starts in 1986.  The scale is chosen so
    that the fitted coefficient is of the same order as the others in the regression.
    """
    v = pd.read_csv(os.path.join(DATA, "FRED", "VIXCLS.csv"), parse_dates=["observation_date"], na_values=".")
    o = pd.read_csv(os.path.join(DATA, "FRED", "VXOCLS.csv"), parse_dates=["observation_date"], na_values=".")
    v = v.set_index("observation_date")["VIXCLS"]; o = o.set_index("observation_date")["VXOCLS"]
    s = v.combine_first(o[o.index < v.index.min()]).dropna()
    q = ((s/100)**2).groupby(s.index.to_period("Q")).mean()
    return 100*q          # percentage points squared / 100, a convenient scale


def main():
    """Run the three specifications over four samples, write the CSV and print the table.

    Where the coefficient on the real rate is positive it is inverted to report an implied Delta, and the
    ratio to the coefficient on the surprise gives an implied gamma - Delta.  For the specification that
    enters the latest growth and the previous expectation separately, the row also carries the test that
    their coefficients sum to zero, which is what aspirations imply.
    """
    d = load()
    d["eps"] = d["g"] - d["mu_lag"]
    d["vix2"] = vix_quarterly().reindex(d.index)
    q = d.index.astype(str)
    zlb = ((q >= "2009Q1") & (q <= "2015Q4")) | ((q >= "2020Q2") & (q <= "2022Q1"))
    samples = {"all": np.ones(len(d), bool),
               "excl. 2020-21": ~(q.str.startswith("2020") | q.str.startswith("2021")),
               "1981-2007": q <= "2007Q4",
               "excl. zero lower bound": ~zlb}
    specs = {"r, eps": ["r", "eps"],
             "r, g, mu_lag": ["r", "g", "mu_lag"],
             "r, eps, vix2": ["r", "eps", "vix2"]}
    rows = []
    for sname, mask in samples.items():
        for spname, xs in specs.items():
            s = d[mask][["mu"] + xs].dropna()
            f = Fit(s, "mu", xs, 4)
            row = dict(sample=sname, spec=spname, n=f.nobs, first=str(s.index[0]), last=str(s.index[-1]), r2=f.rsquared)
            for k in xs:
                row[k] = f.params[k]; row[k + "_se"] = f.bse[k]
            if "eps" in xs and f.params["r"] > 0:
                row["implied_Delta"] = 1/f.params["r"]; row["implied_gamma_minus_Delta"] = f.params["eps"]/f.params["r"]
            if spname == "r, g, mu_lag":
                est, p = f.lin([{"g": 1.0, "mu_lag": 1.0}])   # aspirations: coefficient on mu_lag = -(coefficient on g)
                row["test_g_plus_mulag"] = est[0]; row["p_g_plus_mulag"] = p
            rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ODATA, "spf_pinned_rate_results.csv"), index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
