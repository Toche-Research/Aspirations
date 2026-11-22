"""Least squares with Newey-West standard errors, and the tests built on it."""
import numpy as np
from scipy import stats

def newey_west(X, e, lags):
    """HAC covariance (Bartlett kernel) of OLS coefficients."""
    n = X.shape[0]
    Xe = X * e[:, None]
    S = Xe.T @ Xe
    for l in range(1, lags + 1):
        w = 1 - l / (lags + 1)
        G = Xe[l:].T @ Xe[:-l]
        S += w * (G + G.T)
    XXi = np.linalg.inv(X.T @ X)
    return XXi @ S @ XXi


class Fit:
    """Least squares of y on a constant and xs, with a Newey-West covariance of `lags` lags.

    Exposes nobs, rsquared, params, bse and pvalues keyed by regressor name, the covariance V, and lin()
    for Wald tests of linear restrictions.  Written out here because statsmodels is not a dependency.
    """

    def __init__(self, s, y, xs, lags):
        self.names = ["const"] + xs
        X = np.column_stack([np.ones(len(s))] + [s[x].values for x in xs])
        Y = s[y].values
        b = np.linalg.lstsq(X, Y, rcond=None)[0]
        e = Y - X @ b
        self.V = newey_west(X, e, lags)
        self.b = b
        self.nobs = len(Y)
        self.rsquared = 1 - e @ e / ((Y - Y.mean()) @ (Y - Y.mean()))
        self.params = dict(zip(self.names, b))
        self.bse = dict(zip(self.names, np.sqrt(np.diag(self.V))))
        self.pvalues = {k: 2 * stats.norm.sf(abs(self.params[k] / self.bse[k])) for k in self.names}

    def lin(self, R):
        """Wald test of R b = 0 (rows of R as dicts name->weight); returns (estimates, p-value)."""
        Rm = np.array([[r.get(k, 0.0) for k in self.names] for r in R])
        est = Rm @ self.b
        W = est @ np.linalg.solve(Rm @ self.V @ Rm.T, est)
        return est, stats.chi2.sf(W, len(R))


def ols(d, y, xs, lags=4):
    """Fit y on xs over the rows of d where every variable is observed; return (Fit, sample)."""
    s = d[[y] + xs].dropna()
    return Fit(s, y, xs, lags), s
