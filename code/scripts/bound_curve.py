"""Bound on gamma - Delta from the volatility of ex ante real rates at the aspiration horizon h (the "bound curve").

Model (Proposition 4 with period length h, and its proof, part (iii)): the log price of a bond maturing j periods ahead loads
the latest consumption surprise eps_t once, so the annualized yield at maturity m = j h is
    y_m = const - (gamma - Delta) eps_t / m,     eps_t = growth surprise over the last period (not annualized).
At the aspiration horizon itself (m = h) the annualized one-period rate is r_h = const - (gamma - Delta) e_t, where
e_t = eps_t/h is the ANNUALIZED growth surprise over a period of length h; with i.i.d. growth sd(eps_t) = sigma sqrt(h) and
sd(r_h) = (gamma - Delta) sigma / sqrt(h) (the Timing paragraph of the paper).  If the observed rate is the model rate plus
noise independent of the surprise, then
    (levels)       sd(r_h) >= (gamma - Delta) sd(e_h)            =>  gamma - Delta <= B(h)  = sd(r_h) / sd(e_h),
    (differences)  sd(r_h - r_h,lag) >= sqrt(2) (gamma - Delta) sd(e_h)  =>  gamma - Delta <= B_d(h) = sd(dr_h) / (sqrt 2 sd(e_h)),
since under the model r_h is white noise (i.i.d. growth) and differencing removes slow-moving noise.  The bound is NOT valid
when the noise is correlated with the surprise, e.g. when policy raises the rate after good news (leaning against the wind):
then observed volatility can be below the model's.  Longer maturities m > h bound gamma - Delta by m sd(y_m)/sd(eps_h),
which is looser than the bound at m = h; the curve therefore uses the rate whose maturity equals h.

Data (offline, in data/): h = 1/4: 3-month T-bill (FRED TB3MS, survey month of quarter s) minus the SPF mean forecast of
GDP-deflator inflation over the bill's life (mean of DPGDP2, DPGDP3), 1968Q4-2026Q3; surprise e = first-release annualized
growth of real consumption in quarter s-1 minus the SPF nowcast made in quarter s-1 (DRCONSUM2), 1981Q4-2026Q3, exactly as in
spf_aspirations_test.py.  h = 1: 1-year Treasury (FRED GS1, February) minus the Michigan median expected inflation over the
next 12 months (FRED MICH, February), 1978-2026; surprise = realized growth over Q1-Q4 of the year (mean of the four
first-release annualized quarterly growth rates) minus the forecast of the same four quarters made in the Q1 survey
(mean of DRCONSUM2-5), 1982-2024; variant with the Q4 survey of the previous year (DRCONSUM3-6, rate in November).
Variant with the surprise relative to the one-quarter-ahead forecast made in quarter s-2 (DRCONSUM3(s-2)): the nowcast
is made in the middle of quarter s-1 and so covers only half a period of news (B biased up); the one-quarter-ahead forecast
covers 1.5 periods (B biased down).  Residualized bound: r and e residualized on mu_t and mu_{t-1} (annual: on the
year's forecast), since with time-varying expected growth r = c + Delta mu_t - (gamma-Delta) e_t + v and the levels
inequality needs the constant not to covary with e.  Differences: sd(dr)/sd(de) is reported as well as sd(dr)/(sqrt2 sd(e)),
which presumes a white-noise surprise.
Variants: sigma from the unconditional sd of growth instead of the surprise (the i.i.d. case: real PCE per capita,
FRED A794RX0Q048SBEA, quarterly annualized and annual); ex post rates (bill minus realized deflator inflation, FRED
A191RI1Q225SBEA, next quarter / next four quarters); sub-samples to 2007Q4 and excluding 2020-22.
Reference: the Mehra-Prescott ex post annual bound 1.59 used in the paper (sd 5.67 / sigma 3.57).

Feeds Figure 2 and the 'bound' paragraph of Section 7.  Writes out/data/bound_curve.csv (one row per variant),
out/data/bound_curve_data_quarterly.csv and out/data/bound_curve_data_annual.csv (the assembled series), and
out/figures/fig_bound_curve.pdf.  Three further diagnostics quoted in the paper, the autocorrelation of the
rate and of the surprise and the correlation of the February rate with the previous year's surprise, are
printed to the console only and are not in the CSV files.  The figure path is relative, so run the script
from code/.  See the header of spf_aspirations_test.py for the required data layout."""
import os
import warnings
import numpy as np
import pandas as pd
from aspirations import style  # selects the pgf backend; must precede pyplot
style.setup()
import matplotlib.pyplot as plt
from aspirations import DATA, ODATA
from aspirations.data import load, q_index

warnings.filterwarnings("ignore")
FIG = os.path.join(os.path.dirname(ODATA), "figures")
os.makedirs(FIG, exist_ok=True)
MP_BOUND = 5.67 / 3.57     # Mehra-Prescott (1985) ex post annual bound used in the paper


def fred(name):
    """Read the FRED series `name` from data/FRED/<name>.csv as a dated Series, dropping gaps."""
    s = pd.read_csv(os.path.join(DATA, "FRED", name + ".csv"), parse_dates=["observation_date"], na_values=".")
    return s.set_index("observation_date")[name].dropna()


def quarterly_data():
    """h = 1/4: ex ante real bill rate and SPF nowcast surprise, quarterly (spf_aspirations_test.load conventions)."""
    d = load()
    d["e"] = d["g"] - d["mu_lag"]                       # surprise about quarter s-1, known at survey s
    d["e1"] = d["g"] - d["mu_next"].shift(2)            # the same relative to the one-quarter-ahead forecast of survey s-2
    infl = fred("A191RI1Q225SBEA"); infl.index = infl.index.to_period("Q")
    d["r_expost"] = d["tbill"] - infl.reindex(d.index + 1).values   # bill in quarter s minus realized inflation in s+1
    pce = fred("A794RX0Q048SBEA"); pce.index = pce.index.to_period("Q")
    gq = 400 * np.log(pce).diff()
    d["g_pce"] = gq.reindex(d.index - 1).values          # growth of quarter s-1 (final data)
    return d


def annual_data():
    """h = 1: ex ante one-year real rate (February) and SPF four-quarter surprise, annual."""
    rc = pd.read_excel(os.path.join(DATA, "RCONSUM", "Mean_RCONSUM_Growth.xlsx")); rc.index = q_index(rc.YEAR, rc.QUARTER)
    fr = pd.read_excel(os.path.join(DATA, "RCONSUM", "rcon_first_second_third.xlsx"), sheet_name="DATA", header=4)
    fr = fr[fr["Date"].astype(str).str.match(r"^\d{4}:Q\d$")]
    fr.index = pd.PeriodIndex([pd.Period(x.replace(":", ""), freq="Q") for x in fr["Date"]])
    first = pd.to_numeric(fr["First"], errors="coerce").fillna(pd.to_numeric(fr["Second"], errors="coerce"))
    q1 = rc[rc.QUARTER == 1]; q4 = rc[rc.QUARTER == 4]
    mu_q1 = pd.Series(q1[["DRCONSUM2", "DRCONSUM3", "DRCONSUM4", "DRCONSUM5"]].mean(axis=1).values, index=q1.YEAR.values)
    mu_q4 = pd.Series(q4[["DRCONSUM3", "DRCONSUM4", "DRCONSUM5", "DRCONSUM6"]].mean(axis=1).values, index=q4.YEAR.values + 1)
    g_year = first.groupby(first.index.year).agg(["mean", "count"])
    g_year = g_year[g_year["count"] == 4]["mean"]           # mean of four first-release annualized quarterly growth rates
    gs1, mich = fred("GS1"), fred("MICH")
    rr = (gs1 - mich).dropna()
    r_feb = pd.Series(rr[rr.index.month == 2].values, index=rr[rr.index.month == 2].index.year)
    r_nov = pd.Series(rr[rr.index.month == 11].values, index=rr[rr.index.month == 11].index.year + 1)
    infl = fred("A191RI1Q225SBEA")
    infl_year = infl.groupby(infl.index.year).mean()
    r_feb_expost = pd.Series(gs1[gs1.index.month == 2].values, index=gs1[gs1.index.month == 2].index.year) - infl_year
    pce = fred("A794RX0Q048SBEA")
    pce_year = pce.groupby(pce.index.year).agg(["mean", "count"]); pce_year = pce_year[pce_year["count"] == 4]["mean"]
    g_pce = 100 * np.log(pce_year).diff()
    return pd.DataFrame({"r_feb": r_feb, "r_nov": r_nov, "r_feb_expost": r_feb_expost, "mu_q1": mu_q1, "mu_q4": mu_q4,
                         "g": g_year, "e_q1": g_year - mu_q1, "e_q4": g_year - mu_q4, "g_pce": g_pce})


def residualize(y, X):
    """Residual of y after least squares on a constant and the columns of X."""
    X = np.column_stack([np.ones(len(y))] + [np.asarray(x, float) for x in X])
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def bound(r, e, label, h, primary=False, kind="ex ante, levels", controls=None):
    """Bound on gamma - Delta from one pairing of a rate series r with a surprise series e.

    Returns a row holding the sample, the two standard deviations and those of their first differences,
    the bound in levels sd(r)/sd(e), the bound in differences on the white-noise assumption
    sd(dr)/(sqrt2 sd(e)) and without it sd(dr)/sd(de), and the autocorrelations and correlation that
    bear on whether the independence assumption behind the bound is tenable.  When `controls` is given,
    r and e are also residualized on them and the resulting bound is reported as well.  `primary` marks
    the two headline pairings, which are the filled circles of Figure 2; `h` is the aspiration horizon
    in years and `kind` groups the variants for plotting.
    """
    s = pd.DataFrame({"r": r, "e": e}).dropna()
    sd_r, sd_e, sd_dr, sd_de = s.r.std(), s.e.std(), s.r.diff().std(), s.e.diff().std()
    row = dict(h=h, measure=label, kind=kind, primary=primary, n=len(s), first=str(s.index[0]), last=str(s.index[-1]),
               sd_rate=sd_r, sd_change_rate=sd_dr, sd_surprise=sd_e, sd_change_surprise=sd_de,
               bound_levels=sd_r / sd_e, bound_differences=sd_dr / (np.sqrt(2) * sd_e), bound_differences_de=sd_dr / sd_de,
               autocorr_rate=s.r.autocorr(), autocorr_surprise=s.e.autocorr(), corr_rate_surprise=s[["r", "e"]].corr().iloc[0, 1])
    if controls is not None:                      # residualize r and e on expected growth (time-varying moments)
        c = pd.DataFrame({"r": r, "e": e, **{f"x{i}": x for i, x in enumerate(controls)}}).dropna()
        X = [c[k].values for k in c.columns if k.startswith("x")]
        rr, ee = residualize(c.r.values, X), residualize(c.e.values, X)
        row.update(corr_mu_surprise=np.corrcoef(c.x0, c.e)[0, 1], n_resid=len(c),
                   bound_levels_residualized=rr.std(ddof=1) / ee.std(ddof=1))
    return row


def main():
    """Compute every variant of the bound, write the three CSV files and draw Figure 2."""
    q, a = quarterly_data(), annual_data()
    rows = []
    # h = 1/4
    pre = q.index <= pd.Period("2007Q4", freq="Q"); ex = ~q.index.year.isin([2020, 2021, 2022])
    rows.append(bound(q.r, q.e, "bill - SPF inflation; SPF nowcast surprise", 0.25, True, controls=(q.mu, q.mu_lag)))
    rows.append(bound(q.r[pre], q.e[pre], "same, to 2007Q4", 0.25, controls=(q.mu[pre], q.mu_lag[pre])))
    rows.append(bound(q.r[ex], q.e[ex], "same, excl. 2020-22", 0.25, controls=(q.mu[ex], q.mu_lag[ex])))
    rows.append(bound(q.r, q.e1, "bill - SPF inflation; surprise vs one-quarter-ahead forecast", 0.25, kind="ex ante, levels (1q-ahead)"))
    rows.append(bound(q.r[pre], q.e1[pre], "same, to 2007Q4", 0.25, kind="ex ante, levels (1q-ahead)"))
    rows.append(bound(q.r[ex], q.e1[ex], "same, excl. 2020-22", 0.25, kind="ex ante, levels (1q-ahead)"))
    rows.append(bound(q.r, q.g_pce, "bill - SPF inflation; sd of PCE growth (i.i.d. sigma)", 0.25, kind="ex ante, i.i.d. sigma"))
    rows.append(bound(q.r[ex], q.g_pce[ex], "same, excl. 2020-22", 0.25, kind="ex ante, i.i.d. sigma"))
    rows.append(bound(q.r_expost, q.e, "bill - realized inflation; SPF nowcast surprise", 0.25, kind="ex post"))
    rows.append(bound(q.r_expost, q.g_pce, "bill - realized inflation; sd of PCE growth, 1947-2026", 0.25, kind="ex post, i.i.d. sigma"))
    # h = 1
    pre_a = a.index <= 2007; ex_a = ~a.index.isin([2020, 2021, 2022])
    rows.append(bound(a.r_feb, a.e_q1, "1-year - Michigan (Feb.); SPF four-quarter surprise (Q1 survey)", 1.0, True, controls=(a.mu_q1,)))
    rows.append(bound(a.r_feb[pre_a], a.e_q1[pre_a], "same, to 2007", 1.0))
    rows.append(bound(a.r_feb[ex_a], a.e_q1[ex_a], "same, excl. 2020-22", 1.0))
    rows.append(bound(a.r_nov, a.e_q4, "1-year - Michigan (Nov.); SPF four-quarter surprise (Q4 survey)", 1.0))
    rows.append(bound(a.r_feb, a.g_pce, "1-year - Michigan (Feb.); sd of annual PCE growth (i.i.d. sigma)", 1.0, kind="ex ante, i.i.d. sigma"))
    rows.append(bound(a.r_feb[ex_a], a.g_pce[ex_a], "same, excl. 2020-22", 1.0, kind="ex ante, i.i.d. sigma"))
    rows.append(bound(a.r_feb_expost, a.e_q1, "1-year - realized inflation; SPF four-quarter surprise", 1.0, kind="ex post"))
    rows.append(bound(a.r_feb_expost, a.g_pce, "1-year - realized inflation; sd of annual PCE growth, 1954-2025", 1.0, kind="ex post, i.i.d. sigma"))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ODATA, "bound_curve.csv"), index=False)
    q.to_csv(os.path.join(ODATA, "bound_curve_data_quarterly.csv")); a.to_csv(os.path.join(ODATA, "bound_curve_data_annual.csv"))
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 70)
    print(res.round(3).to_string(index=False))
    qs = q[["r", "e"]].dropna()
    print("bound sample, quarterly: autocorrelation of the ex ante rate", round(qs.r.autocorr(), 3),
          "; of the surprise", round(qs.e.autocorr(), 3), "; corr(rate, latest surprise)", round(qs.corr().iloc[0, 1], 3),
          "; corr(mu, e)", round(q[["mu", "e"]].dropna().corr().iloc[0, 1], 3))
    la = pd.DataFrame({"r": a.r_feb, "e_latest": a.e_q1.shift(1)}).dropna()   # February rate of y vs the surprise of y-1
    print("annual: corr(February rate, surprise of the previous year)", round(la.corr().iloc[0, 1], 3), "n =", len(la),
          "; autocorrelation of the surprise", round(a[["r_feb", "e_q1"]].dropna().e_q1.autocorr(), 3),
          "; corr(mu, e)", round(a[["mu_q1", "e_q1"]].dropna().corr().iloc[0, 1], 3))

    # figure
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    prim = res[res.primary]
    ax.plot(prim.h, prim.bound_levels, "o", color="k", ms=6, label="ex ante real rate, levels")
    ax.plot(prim.h, prim.bound_differences_de, "s", color="k", ms=6, mfc="white", label="ex ante real rate, first differences")
    xp = res[res.kind == "ex post"]
    ax.plot(xp.h, xp.bound_levels, "x", color="k", ms=7, mew=1.2, zorder=3, label="ex post real rate, levels")
    ax.axhline(MP_BOUND, color="k", ls="--", lw=1, label="annual bill returns, 1889\u20131978 (Table 1)")
    sub = res[(res.kind == "ex ante, levels") & ~res.primary]
    ax.plot(sub.h, sub.bound_levels, marker="_", color="k", ms=9, mew=1, ls="none", label="levels: sub-samples, other timing")
    ax.set_xlim(0, 1.1); ax.set_ylim(0, 3.2)
    ax.set_xticks([0.25, 0.5, 0.75, 1.0]); ax.set_xticklabels(["1/4", "1/2", "3/4", "1"])
    ax.set_xlabel("aspiration horizon $h$ (years)"); ax.set_ylabel(r"bound on $\gamma-\Delta$")
    # the fractional sizes here are squeezed to fit this figure; style.ANNOT is the default
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_bound_curve.pdf"), bbox_inches="tight")
    return res


if __name__ == "__main__":
    main()
