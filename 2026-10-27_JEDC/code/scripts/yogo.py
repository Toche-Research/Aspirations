"""Replication of Yogo (2004, REStat 86(3), Table: EIS estimates with the risk-free rate) and the aspirations test.

Data: Campbell (2003) international data as distributed by Yogo, in data/yogo2004/<CTRY>Q.txt (quarterly) and
<CTRY>Y.txt (annual). Columns: DATE r dp rf inf dc rr rrf z1..z4, with z1..z4 the instruments (dp, rf, inf, dc) lagged
twice. As in Yogo's ivreg.prg: drop rows with missing values; for the short US quarterly sample start at 1970.3. Estimators (k-class with the constant projected out, homoskedastic standard errors with T-2 degrees of
freedom), exactly as in Yogo's tsls.prg, fuller.prg and liml.prg.

The data files are not stored unpacked in the repository: before running, create the data directory and
unpack Yogo's archive into it, from the repository root,
    the repository ships these files already, unzipped from Yogo's EIS_Data.zip

Feeds no numbered table or figure.  It establishes that the replication reproduces Yogo's published
estimates, so that the aspirations test run on the same data in yogo_aspirations_test.py rests on a
verified reading of them.  All results are printed; out/yogo_replication.txt is the saved console output,
    python3 yogo.py > out/yogo_replication.txt"""
from pathlib import Path
import numpy as np

DATA = Path(__file__).resolve().parents[1] / 'data' / 'yogo2004'
SERIES = [('USA', 'Q', False), ('AUL', 'Q', False), ('CAN', 'Q', False), ('FR', 'Q', False), ('GER', 'Q', False),
          ('ITA', 'Q', False), ('JAP', 'Q', False), ('NTH', 'Q', False), ('SWD', 'Q', False), ('SWT', 'Q', False),
          ('UK', 'Q', False), ('USA', 'Q', True), ('SWD', 'Y', False), ('UK', 'Y', False), ('USA', 'Y', False)]


def load(ctry, freq, short=False):
    """Read one country-frequency data file and return its columns as a dict of arrays.

    Missing values are written '.' in Yogo's files; rows with any missing value are dropped, as in his
    ivreg.prg.  With short=True the sample starts at the 95th data row, which is 1970.3 for the United
    States quarterly series.
    """
    rows = []
    for line in (DATA / f'{ctry}{freq}.txt').read_text().splitlines()[1:]:
        f = line.split()
        if not f:
            continue
        rows.append([np.nan if v == '.' else float(v) for v in f])
    d = np.array(rows)
    # Gauss reads the header line as a row of missing values, so its row 96 is data row 95 (0-based 94), and its
    # data[2:rows] keeps every data row; rows with missing instruments are then dropped.
    d = d[94:] if short else d
    d = d[~np.isnan(d).any(axis=1)]
    cols = 'DATE r dp rf inf dc rr rrf z1 z2 z3 z4'.split()
    return {c: d[:, i] for i, c in enumerate(cols)}


def demean(a):
    """Subtract column means, treating a one-dimensional input as a single column.

    Demeaning is how the constant is projected out of the k-class estimators below.
    """
    a = np.atleast_2d(a.T).T if a.ndim == 1 else a
    return a - a.mean(axis=0)


def kclass(y, X, Z, kind='tsls'):
    """k-class estimator of y on X (n x p) with instruments Z (n x k2); constant projected out. Returns beta, se."""
    y = demean(y)[:, 0]; X = demean(X); Z = demean(Z)
    T, k2 = Z.shape
    PZ = Z @ np.linalg.solve(Z.T @ Z, Z.T)
    MZ = np.eye(T) - PZ
    if kind == 'tsls':
        k = 1.0
    else:
        Yb = np.column_stack([y, X])
        kl = np.min(np.real(np.linalg.eigvals(np.linalg.solve(Yb.T @ MZ @ Yb, Yb.T @ Yb))))
        k = kl if kind == 'liml' else kl - 1.0 / (T - k2 - 1)
    W = np.eye(T) - k * MZ
    A = X.T @ W @ X
    beta = np.linalg.solve(A, X.T @ W @ y)
    u = y - X @ beta
    varu = u @ u / (T - 1 - X.shape[1])
    se = np.sqrt(np.diag(varu * np.linalg.inv(A)))
    return beta, se


def first_stage_F(x, Z):
    """First-stage F statistic of the regression of x on the instruments Z, constant projected out.

    The conventional weak-instrument diagnostic: a low value means Z explains little of x.
    """
    x = demean(x)[:, 0]; Z = demean(Z); T, k2 = Z.shape
    b = np.linalg.lstsq(Z, x, rcond=None)[0]; e = x - Z @ b
    r2 = 1 - e @ e / (x @ x)
    return (r2 / k2) / ((1 - r2) / (T - k2 - 1))


if __name__ == '__main__':
    print('Replication of Yogo (2004), risk-free rate: 1/psi from rrf on dc; psi from dc on rrf (TSLS, Fuller, LIML)')
    for ctry, freq, short in SERIES:
        d = load(ctry, freq, short)
        Z = np.column_stack([d['z1'], d['z2'], d['z3'], d['z4']])
        out = []
        for kind in ('tsls', 'fuller', 'liml'):
            b, s = kclass(d['rrf'], d['dc'][:, None], Z, kind); out.append(f'{b[0]:7.2f} ({s[0]:5.2f})')
        for kind in ('tsls', 'fuller', 'liml'):
            b, s = kclass(d['dc'], d['rrf'][:, None], Z, kind); out.append(f'{b[0]:6.2f} ({s[0]:4.2f})')
        print(f'{ctry:4s}{freq} {d["DATE"][0]:.1f}-{d["DATE"][-1]:.1f} T={len(d["DATE"]):3d}  ' + '  '.join(out))
