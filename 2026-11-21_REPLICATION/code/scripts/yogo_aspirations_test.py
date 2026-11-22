"""Test of aspirations against catching up (Proposition 'Time-varying moments'), on Yogo's (2004) data.

Proposition 5(ii)-(iii): with aspirations, r^f_{t+1} = c + Delta mu_t + (gamma - Delta)(mu_{t-1} - g_t) + variance terms,
where mu_{t-1} - g_t = -eps_t is last period's growth surprise; with catching up, r^f_{t+1} = c + gamma mu_t - (gamma-Delta) g_t.
In the data (Yogo's timing: rrf_{t+1} and dc_{t+1} over period t+1), regress rrf_{t+1} on dc_{t+1} and dc_t with
instruments dated t-1 (Yogo's z: dp, rf, inf, dc lagged twice):
    aspirations:  coefficient on dc_t = 0,              on dc_{t+1} = Delta
    catching up:  coefficient on dc_t = -(gamma-Delta) < 0, on dc_{t+1} = gamma
    power:        coefficient on dc_t = 0,              on dc_{t+1} = gamma
Robustness: instruments lagged once more (dated t-2), against time aggregation of consumption.
Inference: LIML and TSLS with Yogo's homoskedastic standard errors, first-stage F statistics, and an Anderson-Rubin
confidence set for the coefficient on dc_t obtained by projection (robust to weak instruments).

Feeds the Yogo-data sentences of Section 7 ('Euler-equation tests'): the first-stage F statistic for
lagged growth, and the Anderson-Rubin sets for its coefficient, which reach the edge of the [-20, 20]
grid in every series at both instrument dates.  Requires the unpacked Yogo data; see yogo.py for the
one command that unpacks them.  All results are printed; out/yogo_aspirations_test.txt is the saved
console output,
    python3 yogo_aspirations_test.py > out/yogo_aspirations_test.txt"""
import numpy as np
from scipy import stats
from yogo import DATA, SERIES, kclass, first_stage_F, demean

COLS = 'DATE r dp rf inf dc rr rrf z1 z2 z3 z4'.split()


def load_raw(ctry, freq):
    """Read one country-frequency data file as a plain array, keeping rows with missing values.

    Unlike yogo.load, nothing is dropped here: the lagged growth and the extra instrument lag are built
    first, and only then are incomplete rows removed, so that no usable observation is lost.
    """
    rows = []
    for line in (DATA / f'{ctry}{freq}.txt').read_text().splitlines()[1:]:
        f = line.split()
        if f:
            rows.append([np.nan if v == '.' else float(v) for v in f])
    return np.array(rows)


def build(ctry, freq, short, extra_lag=False):
    """Assemble the regression sample for one series.

    Returns (dates, rrf, [dc_{t+1}, dc_t], Z).  Yogo's timing has rrf and dc dated t+1 on each row, so
    lagged growth dc_t is the previous row's dc.  Z is Yogo's z1..z4, dated t-1, lagged once more when
    extra_lag is set, which guards against time aggregation of consumption.  Rows with any missing value
    are dropped last, after the lags are formed.
    """
    d = load_raw(ctry, freq)
    i = {c: k for k, c in enumerate(COLS)}
    lag = lambda col, n: np.concatenate([np.full(n, np.nan), d[:-n, i[col]]])
    dc_lag = lag('dc', 1)
    Z = d[:, [i['z1'], i['z2'], i['z3'], i['z4']]]
    if extra_lag:
        Z = np.vstack([np.full((1, 4), np.nan), Z[:-1]])
    M = np.column_stack([d[:, i['DATE']], d[:, i['rrf']], d[:, i['dc']], dc_lag, Z])
    if short:
        M = M[94:]
    M = M[~np.isnan(M).any(axis=1)]
    return M[:, 0], M[:, 1], M[:, 2:4], M[:, 4:]


def ar_projection_ci(y, X, Z, grid2, grid1, alpha=0.05):
    """Anderson-Rubin set for (b1, b2), projected on b2: b2 is in the set if some b1 is not rejected."""
    yz, Xz, Zz = demean(y)[:, 0], demean(X), demean(Z)
    T, k2 = Zz.shape
    crit = stats.f.ppf(1 - alpha, k2, T - k2 - 1)
    PZ = Zz @ np.linalg.solve(Zz.T @ Zz, Zz.T)
    keep = []
    for b2 in grid2:
        for b1 in grid1:
            e = yz - Xz @ np.array([b1, b2])
            num = e @ PZ @ e / k2
            den = (e @ e - e @ PZ @ e) / (T - k2 - 1)
            if num / den < crit:
                keep.append(b2)
                break
    if not keep:
        return 'empty'
    keep = np.array(keep)
    lo, hi = keep.min(), keep.max()
    edge = (lo <= grid2[0] + 1e-12) or (hi >= grid2[-1] - 1e-12)
    return f'[{lo:6.2f},{hi:6.2f}]' + ('*' if edge else ' ')


if __name__ == '__main__':
    grid2 = np.linspace(-20, 20, 801); grid1 = np.linspace(-40, 40, 801)
    for extra in (False, True):
        print('Instruments dated t-1 (Yogo)' if not extra else '\nInstruments dated t-2 (robust to time aggregation)')
        print('series        T   F(dc_t+1) F(dc_t)  | LIML: dc_t+1        dc_t          | TSLS dc_t       | AR 95% set for dc_t coef')
        for ctry, freq, short in SERIES:
            dates, y, X, Z = build(ctry, freq, short, extra)
            F1, F2 = first_stage_F(X[:, 0], Z), first_stage_F(X[:, 1], Z)
            bl, sl = kclass(y, X, Z, 'liml'); bt, st = kclass(y, X, Z, 'tsls')
            ci = ar_projection_ci(y, X, Z, grid2, grid1)
            tag = f'{ctry}{freq}{" short" if short else ""}'
            print(f'{tag:11s} {len(y):4d}  {F1:7.2f} {F2:8.2f}  | {bl[0]:7.2f} ({sl[0]:6.2f}) {bl[1]:7.2f} ({sl[1]:6.2f}) | '
                  f'{bt[1]:6.2f} ({st[1]:5.2f}) | {ci}')
    print('\n* : the set reaches the edge of the grid [-20, 20] (unbounded in practice).')
