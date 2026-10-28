"""Numerical checks of the four derivations of Appendix B of 'Ex Ante Aspirations'.

(1) Individual Euler equation with aspirations set one period ahead and a riskless return R = e^r (isoelastic example,
    log c_{t+1} = log c_t + mu_t + eps_{t+1}, eps ~ N(0, s2)):  Delta mu_t = log beta + r + gamma^2 s2 / 2 + (gamma - Delta) eps_t,
    where eps_t = log c_t - log c_{t-1} - mu_{t-1} is the current surprise.  Hence consumption growth
    g_{t+1} = mu_t + eps_{t+1} = kappa_t + theta_h eps_t + eps_{t+1} with theta_h = (gamma - Delta)/Delta: MA(1).
    Checked by solving the exact Euler equation U_c(c_t, a_t) = beta R E_t U_c(c_{t+1}, a_{t+1}) for mu_t by bisection,
    with the expectation by Gauss-Hermite quadrature and the aspirations from the first-order condition
    E U_a = -B'(a), for random parameters.
(2) Satisfaction S = -log(1/(gamma-1) - U(c_{t+1}, a_{t+1})) = const + (Delta-1)(log c_t + mu_t) + (gamma-1) eps_{t+1}
    + (gamma-Delta)(gamma-1) s2 / 2, checked by finite differences.
(3) Under certainty U(c, a*(c)) = 1/(gamma-1) - K c^{1-Delta}: falls with c when Delta < 1.
(4) What an IV regression of annual growth on lagged annual growth identifies under the MA(1), by
    simulation, including the triangular instrument loadings implied by Dynan's (2000) design; see the
    long comment before alpha_formula for the derivation.

Feeds the 'individual Euler equation and habits' paragraph of Section 7 and Appendix B: the loadings omega
and omega_- for January and June interviews, the resulting upper limits on theta_h and on gamma - Delta,
and the finding that a quarterly aspiration horizon leaves the design without power.  Checks (1) to (3)
are printed; check (4) is printed and also written to out/micro_euler_dynan.txt.
Run: python3 checks_micro_euler.py"""
import numpy as np

rng = np.random.default_rng(0)
x, w = np.polynomial.hermite.hermgauss(60)


def E_normal(f, m, s):
    """E f(z), z ~ N(m, s^2), by Gauss-Hermite."""
    return (w * f(m + np.sqrt(2) * s * x)).sum() / np.sqrt(np.pi)


def aspiration(logc, mu, s2, gamma, theta, b, kappa):
    """(theta+b) log a = log(kappa(gamma-1)/(1+theta)) - log E c^{1-gamma}, E over log c ~ N(logc+mu, s2)."""
    logE = (1 - gamma) * (logc + mu) + (1 - gamma) ** 2 * s2 / 2
    return (np.log(kappa * (gamma - 1) / (1 + theta)) - logE) / (theta + b)


def check_euler(trials=20):
    """Check the individual Euler equation at random parameters; return the largest error in mu_t.

    Expected growth is solved from the exact Euler equation by bisection, with the expectation taken by
    Gauss-Hermite quadrature and the aspiration from its first-order condition, and compared with the
    closed form Delta mu_t = log beta + r + gamma^2 sigma^2/2 + (gamma - Delta) eps_t.
    """
    worst = 0.0
    for _ in range(trials):
        gamma = rng.uniform(1.2, 4); theta = rng.uniform(gamma - 1 + 0.2, 8); b = rng.uniform(0, 2); kappa = rng.uniform(0.5, 3)
        beta = rng.uniform(0.9, 0.99); r = rng.uniform(-0.02, 0.08); s = rng.uniform(0.005, 0.1); s2 = s * s
        Delta = 1 + (gamma - 1) * (b - 1) / (theta + b)
        logc_1, mu_1, eps = rng.normal(0, 1), rng.uniform(-0.05, 0.05), rng.normal(0, s)
        logc = logc_1 + mu_1 + eps                                 # log c_t
        loga = aspiration(logc_1, mu_1, s2, gamma, theta, b, kappa)   # a_t, set at t-1
        Uc_t = np.exp(-gamma * logc + (1 + theta) * loga)

        def rhs(mu):
            loga1 = aspiration(logc, mu, s2, gamma, theta, b, kappa)
            return beta * np.exp(r) * E_normal(lambda lc1: np.exp(-gamma * lc1 + (1 + theta) * loga1), logc + mu, s)

        lo, hi = -2.0, 2.0
        for _ in range(200):                                       # rhs is decreasing in mu
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if rhs(mid) > Uc_t else (lo, mid)
        mu_num = (lo + hi) / 2
        mu_formula = (np.log(beta) + r + gamma ** 2 * s2 / 2 + (gamma - Delta) * eps) / Delta
        worst = max(worst, abs(mu_num - mu_formula))
    return worst


def check_satisfaction(trials=20):
    """Check the satisfaction measure by finite differences; return the largest error.

    Confirms that S responds to anticipated log consumption and to expected growth with coefficient
    Delta - 1, to the surprise with gamma - 1, and to the conditional variance with (gamma-Delta)(gamma-1)/2.
    """
    worst = 0.0
    for _ in range(trials):
        gamma = rng.uniform(1.2, 4); theta = rng.uniform(gamma - 1 + 0.2, 8); b = rng.uniform(0, 2); kappa = rng.uniform(0.5, 3)
        s2 = rng.uniform(0.005, 0.1) ** 2; Delta = 1 + (gamma - 1) * (b - 1) / (theta + b)
        logc, mu, eps1 = rng.normal(), rng.uniform(-0.05, 0.05), rng.normal(0, 0.05)

        def S(logc, mu, eps1, s2):
            loga1 = aspiration(logc, mu, s2, gamma, theta, b, kappa)
            U = (np.exp((1 - gamma) * (logc + mu + eps1) + (1 + theta) * loga1) - 1) / (1 - gamma)
            return -np.log(1 / (gamma - 1) - U)

        h = 1e-6
        d_anticipated = (S(logc + h, mu, eps1, s2) - S(logc - h, mu, eps1, s2)) / (2 * h)
        d_mu = (S(logc, mu + h, eps1, s2) - S(logc, mu - h, eps1, s2)) / (2 * h)
        d_surprise = (S(logc, mu, eps1 + h, s2) - S(logc, mu, eps1 - h, s2)) / (2 * h)
        d_s2 = (S(logc, mu, eps1, s2 + h) - S(logc, mu, eps1, s2 - h)) / (2 * h)
        worst = max(worst, abs(d_anticipated - (Delta - 1)), abs(d_mu - (Delta - 1)), abs(d_surprise - (gamma - 1)),
                    abs(d_s2 - (gamma - Delta) * (gamma - 1) / 2))
    return worst


def check_certainty(trials=20):
    """Check that, under certainty, log(1/(gamma-1) - U(c, a*(c))) has elasticity 1 - Delta in c."""
    worst = 0.0
    for _ in range(trials):
        gamma = rng.uniform(1.2, 4); theta = rng.uniform(gamma - 1 + 0.2, 8); b = rng.uniform(0, 2); kappa = rng.uniform(0.5, 3)
        Delta = 1 + (gamma - 1) * (b - 1) / (theta + b)
        c = np.exp(rng.normal())
        a = ((kappa * (gamma - 1) / (1 + theta)) * c ** (gamma - 1)) ** (1 / (theta + b))    # U_a + B' = 0
        U = lambda c, a: (c ** (1 - gamma) * a ** (1 + theta) - 1) / (1 - gamma)
        # d log(1/(gamma-1) - U(c, a*(c))) / d log c should be 1 - Delta
        h = 1e-6
        f = lambda lc: np.log(1 / (gamma - 1) - U(np.exp(lc), ((kappa * (gamma - 1) / (1 + theta)) * np.exp(lc) ** (gamma - 1)) ** (1 / (theta + b))))
        worst = max(worst, abs((f(np.log(c) + h) - f(np.log(c) - h)) / (2 * h) - (1 - Delta)))
    return worst




# ----------------------------------------------------------------------------------------------------------------------
# (4) What an IV regression of annual growth on lagged annual growth identifies under the MA(1)
#     when the instrument is correlated with the surprises of the regressor's own year (Dynan 2000, AER 90(3), 391-406).
#
#     Sub-period (aspiration-horizon) growth: g_s = kappa + theta_h eps_{s-1} + eps_s, n sub-periods per year.
#     Annual growth over year t (sub-periods n(t-1)+1, ..., nt):
#         G_t = n kappa + theta_h l_{t-1} + (1 + theta_h) F_t + l_t,
#     with l_t the surprise of the last sub-period of year t and F_t the sum of the first n-1.  For an instrument z
#     orthogonal to the surprises of year t+1, the IV estimand of G_{t+1} on G_t is
#         alpha = theta_h Cov(z, l_t) / [theta_h Cov(z, l_{t-1}) + (1 + theta_h) Cov(z, F_t) + Cov(z, l_t)]
#               = theta_h w / (1 - w_ + theta_h (1 - w)),
#     where w = Cov(z, l_t)/S, w_ = Cov(z, l_{t-1})/S, S = Cov(z, l_{t-1}) + Cov(z, F_t) + Cov(z, l_t).
#     Equal loading on the n sub-periods of year t and none on year t-1: alpha = theta_h / (n + (n-1) theta_h).
#     Checked below by simulation with a linear instrument and with Dynan-type range dummies.  Dynan's growth
#     instruments (growth of calendar-year income and hours totals) load the months of two calendar years with
#     triangular weights under random-walk income: weights_triangle computes (w, w_), simulate_iv_triangle checks them.


def alpha_formula(theta_h, w, w_):
    """IV estimand of lagged on current annual growth, equation (12) of the paper.

    w is the share of the instrument's covariance with the surprises that falls in the last horizon of the
    lagged year, w_ the share in the horizon just before that year.
    """
    return theta_h * w / (1 - w_ + theta_h * (1 - w))


def simulate_iv(theta_h, n, weights, years=60, households=20000, noise_sd=1.0, dummies=False, seed=1):
    """IV coefficient of G_{t+1} on G_t, year effects removed, instrument z_t = sum_k weights[k] eps_{(k-th sub-period)}
    + noise, where weights has n+1 entries: index 0 is the last sub-period of year t-1, indices 1..n the sub-periods
    of year t.  Sub-period surprises i.i.d. N(0,1).  Time effects absorbed by demeaning within year (kappa = 0 here)."""
    r = np.random.default_rng(seed)
    S = n * years
    eps = r.standard_normal((households, S + 1))                       # eps[:, s], s = 0..S; sub-period s of the panel
    g = theta_h * eps[:, :-1] + eps[:, 1:]                             # g_s for s = 1..S
    G = g.reshape(households, years, n).sum(axis=2)                    # annual growth, years 0..years-1
    # instrument for regressor G_t (year t, 0-indexed): eps columns n t .. n t + n (column n t = last sub-period of year t-1)
    z = np.zeros((households, years))
    for t in range(1, years):
        block = eps[:, n * t: n * t + n + 1]                            # n+1 columns: last of t-1, then year t
        z[:, t] = block @ np.asarray(weights) + noise_sd * r.standard_normal(households)
    y = G[:, 2:].ravel(); x = G[:, 1:-1].ravel(); zz = z[:, 1:-1].ravel()
    # year effects: demean within year (balanced panel)
    def demean(v):
        v = v.reshape(households, years - 2); return (v - v.mean(axis=0)).ravel()
    y, x, zz = demean(y), demean(x), demean(zz)
    if dummies:                                                        # Dynan-type: dummies for ranges of z
        qs = np.quantile(zz, [1 / 8, 2 / 8, 3 / 8, 4 / 8, 5 / 8, 6 / 8, 7 / 8])
        Z = np.column_stack([(zz > lo) & (zz <= hi) for lo, hi in zip(np.r_[-np.inf, qs], np.r_[qs, np.inf])]).astype(float)[:, 1:]
        Z = Z - Z.mean(axis=0)
        xhat = Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
        return (xhat @ y) / (xhat @ x)
    return (zz @ y) / (zz @ x)


def weights_triangle(M, n):
    """Shares (w, w_) of eq. (12) for Dynan's growth instruments, growth of calendar-year income (or hours) totals, when
    income is a monthly random walk whose innovations are proportional to the consumption surprises.  Consumption is
    measured in interview month M (1..12), so the lagged period runs from month M of calendar year A to month M of A+1,
    and the instrument is Y_A - Y_{A-1}: a total over year A loads an innovation in its month m with (13-m)/12 and one in
    month m of A-1 with (m-1)/12, so the growth loads month m of A with (13-m)/12, month m of A-1 with (m-1)/12, and
    nothing in A+1 (triangular weights).  n horizons per year of 12/n months; the last horizon of the lagged period is
    its last 12/n months, the previous period runs from month M of A-1 to month M of A."""
    hor = 12 // n

    def wgt(year_off, m):                 # year_off: -1 for A-1, 0 for A, +1 for A+1
        return (m - 1) / 12 if year_off == -1 else (13 - m) / 12 if year_off == 0 else 0.0

    lag = [wgt(0, M + k) if M + k <= 12 else wgt(1, M + k - 12) for k in range(12)]        # lagged period, month by month
    prev = [wgt(-1, M + k) if M + k <= 12 else wgt(0, M + k - 12) for k in range(12)]      # period before it
    l_t, F_t, l_prev = sum(lag[12 - hor:]), sum(lag[:12 - hor]), sum(prev[12 - hor:])
    S = l_t + F_t + l_prev
    return l_t / S, l_prev / S


def simulate_iv_triangle(theta_h, n, M, years=80, households=20000, seed=3):
    """IV coefficient with monthly random-walk income, consumption surprises equal to the monthly income innovations
    summed over each horizon, consumption measured in month M, instrument = growth of the calendar-year average of
    income over the two calendar years ending before the dependent variable's period (Dynan's dating)."""
    r = np.random.default_rng(seed)
    hor = 12 // n
    u = r.standard_normal((households, years * 12))                   # monthly innovations, calendar months 1..12
    Y = np.cumsum(u, axis=1).reshape(households, years, 12).mean(axis=2)   # calendar-year average income
    P = years - 2                                                     # period p starts in month M of calendar year p
    e = np.stack([u[:, p * 12 + M - 1: p * 12 + M + 11].reshape(households, n, hor).sum(axis=2) for p in range(P)], axis=1)
    e = e.reshape(households, P * n)                                  # horizon surprises in order
    g = theta_h * np.concatenate([np.zeros((households, 1)), e[:, :-1]], axis=1) + e
    G = g.reshape(households, P, n).sum(axis=2)
    z = np.zeros((households, P)); z[:, 1:] = Y[:, 1:P] - Y[:, 0:P - 1]   # for regressor G_p: Y_p - Y_{p-1}
    dm = lambda v: (v - v.mean(axis=0)).ravel()
    y, x, zz = dm(G[:, 2:]), dm(G[:, 1:-1]), dm(z[:, 1:-1])
    return (zz @ y) / (zz @ x)


def theta_bound(a_up, w, w_):
    """Upper limit on theta_h from an upper limit a_up on alpha; none if a_up >= w/(1-w)."""
    return a_up * (1 - w_) / (w - a_up * (1 - w)) if a_up < w / (1 - w) else float("inf")


def check_dynan_iv(out_path=None):
    """Check the IV estimand against simulation and report the bounds it places on theta_h.

    First a set of illustrative instrument loadings, each checked twice, with the instrument entered
    linearly and as Dynan-style range dummies.  Then the triangular loadings that growth of calendar-year
    income totals implies under random-walk income, for annual and quarterly aspiration horizons and
    interview months January to June, each with the upper limit on theta_h and on gamma - Delta that
    Dynan's estimate allows.  Writes the text to out_path when given, and returns the largest discrepancy
    between simulation and formula.
    """
    lines = []
    cases = [  # illustrative loadings: (label, theta_h, n, weights over [last of t-1, sub-periods 1..n of year t])
        ("illustrative: annual aspirations, instrument loads the lagged period only", 0.3, 1, [0, 1]),
        ("illustrative: annual, loads 7/12 on the lagged period and 5/12 on the one before", 0.3, 1, [5 / 12, 7 / 12]),
        ("illustrative: quarterly aspirations, equal loading on the four quarters, none earlier", 0.3, 4, [0, 1, 1, 1, 1]),
        ("illustrative: quarterly, loads first three quarters and the last of the previous year", 0.3, 4, [1, 1, 1, 1, 0]),
        ("illustrative: quarterly, unequal loading on the four quarters", 0.3, 4, [0, 0.5, 1, 1, 0.25]),
        ("illustrative: monthly aspirations, equal loading", 0.3, 12, [0] + [1] * 12),
        ("illustrative: quarterly, equal loading, large theta_h", 1.0, 4, [0, 1, 1, 1, 1]),
    ]
    worst = 0.0
    for label, th, n, wts in cases:
        Ssum = sum(wts); w = wts[-1] / Ssum; w_ = wts[0] / Ssum
        formula = alpha_formula(th, w, w_)
        sim = simulate_iv(th, n, wts)
        simd = simulate_iv(th, n, wts, dummies=True, seed=2)
        worst = max(worst, abs(sim - formula), abs(simd - formula))
        lines.append(f"{label}: theta_h={th}, n={n}, w={w:.4f}, w_={w_:.4f}: formula {formula:.4f}, "
                     f"simulated IV {sim:.4f}, with range dummies {simd:.4f}")
    # Dynan's instruments: growth of calendar-year totals, triangular loading under random-walk income
    a_up = -0.039 + 1.96 * 0.069
    lines.append(f"Dynan (2000) Table 2 col. 1: alpha = -0.039 (0.069); upper 95% limit {a_up:.4f}")
    lines.append("Triangular loading (growth of calendar-year totals, random-walk income), interview month M:")
    worst_tri = 0.0
    for n in (1, 4):
        for M in (1, 2, 3, 4, 5, 6):
            w, w_ = weights_triangle(M, n)
            th = 0.3
            f = alpha_formula(th, w, w_); sim = simulate_iv_triangle(th, n, M)
            worst_tri = max(worst_tri, abs(sim - f))
            b = theta_bound(a_up, w, w_)
            bstr = (f"theta_h <= {b:.3f}; gamma-Delta <= {0.5 * b:.3f} (Delta=1/2), {0.9 * b:.3f} (Delta=0.9)"
                    if np.isfinite(b) else f"no bound (max alpha = w/(1-w) = {w / (1 - w):.3f} < {a_up:.3f})")
            lines.append(f"  n={n}, M={M}: w={w:.3f}, w_={w_:.3f}; formula (theta_h=0.3) {f:.4f}, simulated {sim:.4f}; {bstr}")
    lines.append(f"max |simulated - formula|: illustrative cases {worst:.4f}, triangular cases {worst_tri:.4f}")
    text = "\n".join(lines)
    print(text)
    if out_path:
        with open(out_path, "w") as f:
            f.write(text + "\n")
    return max(worst, worst_tri)

if __name__ == "__main__":
    import os
    print("Euler equation, max |mu numerical - formula|:", check_euler())
    print("satisfaction derivatives, max error:", check_satisfaction())
    print("U(c,a*(c)) exponent 1-Delta, max error:", check_certainty())
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out", "micro_euler_dynan.txt")
    check_dynan_iv(out)
