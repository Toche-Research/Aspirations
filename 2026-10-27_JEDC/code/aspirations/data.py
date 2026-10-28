"""Assembly of the quarterly data set from the survey workbooks and the FRED series."""
import os

import numpy as np
import pandas as pd

from . import DATA

def q_index(year, quarter):
    """Build a quarterly PeriodIndex from parallel arrays of years and quarter numbers."""
    return pd.PeriodIndex([pd.Period(year=int(y), quarter=int(q), freq="Q") for y, q in zip(year, quarter)])


def load():
    """Assemble the quarterly regression sample and return it as a DataFrame indexed by survey quarter.

    Columns: mu and mu_next, the survey's nowcast and one-quarter-ahead forecast of real consumption
    growth; their one-survey lags mu_lag and mu_next_lag; pi_e, expected GDP-deflator inflation over the
    bill's life; tbill, the three-month rate in the survey month; g, first-release growth of the previous
    quarter, which panellists knew when they answered, with g_filled marking the quarters where the second
    release stands in for a missing first release; disp and disp_lag, the cross-sectional dispersion of the
    forecast; and r = tbill - pi_e, the ex ante real rate.  All variables are in annualized percentage points.
    """
    rc = pd.read_excel(os.path.join(DATA, "RCONSUM", "Mean_RCONSUM_Growth.xlsx"))
    rc.index = q_index(rc.YEAR, rc.QUARTER)
    pg = pd.read_excel(os.path.join(DATA, "RCONSUM", "PGDP", "Mean_PGDP_Growth.xlsx"))
    pg.index = q_index(pg.YEAR, pg.QUARTER)

    fr = pd.read_excel(os.path.join(DATA, "RCONSUM", "rcon_first_second_third.xlsx"), sheet_name="DATA", header=4)
    fr = fr[fr["Date"].astype(str).str.match(r"^\d{4}:Q\d$")]
    fr.index = pd.PeriodIndex([pd.Period(d.replace(":", ""), freq="Q") for d in fr["Date"]])
    first = pd.to_numeric(fr["First"], errors="coerce")
    second = pd.to_numeric(fr["Second"], errors="coerce")
    g_rel = first.fillna(second)          # first release; second release where the first is missing
    g_fill = first.isna() & second.notna()

    disp = pd.read_excel(os.path.join(DATA, "RCONSUM", "Dispersion_RCONSUM.xlsx"), sheet_name="D2", header=None)
    disp = disp.iloc[10:, [0, 3]].dropna(subset=[0])
    disp.index = pd.PeriodIndex([pd.Period(str(s), freq="Q") for s in disp[0]])
    d2 = pd.to_numeric(disp[3], errors="coerce")

    tb = pd.read_csv(os.path.join(DATA, "FRED", "TB3MS.csv"), parse_dates=["observation_date"])
    tb = tb[tb.observation_date.dt.month.isin([2, 5, 8, 11])]   # middle month of the quarter (survey month)
    tb.index = tb.observation_date.dt.to_period("Q")
    tb = tb["TB3MS"]

    d = pd.DataFrame({
        "mu": rc["DRCONSUM2"],
        "mu_next": rc["DRCONSUM3"],
        "pi_e": (pg["DPGDP2"] + pg["DPGDP3"]) / 2,
        "tbill": tb,
    })
    d["mu_lag"] = d["mu"].shift(1)
    d["mu_next_lag"] = d["mu_next"].shift(1)
    d["g"] = g_rel.reindex(d.index - 1).values          # quarter s-1, released before survey s
    d["g_filled"] = g_fill.reindex(d.index - 1).values
    d["disp"] = d2.reindex(d.index)
    d["disp_lag"] = d["disp"].shift(1)
    d["r"] = d["tbill"] - d["pi_e"]
    return d
