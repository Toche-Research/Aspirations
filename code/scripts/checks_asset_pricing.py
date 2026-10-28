"""Checks of the asset-pricing results of 'Ex Ante Aspirations'.

Part A (sympy): the long-run elasticity of substitution of an anticipation-stock model with homogeneous utility
                (Gomez and Monteiro, Macroeconomic Dynamics, Proposition 1) is the inverse of 1 - (degree of homogeneity).
Part B (numpy/scipy): aspirations set one period ahead, isoelastic example, lognormal consumption growth.
    B1  the aspiration first-order condition E_t U_a(c', a') = -kappa against the closed form;
    B2  the stochastic discount factor and the risk-free rate against direct quadrature;
    B3  the term structure of strip risk premia against Monte Carlo;
    B4  the price-consumption ratio of the consumption claim against its recursion.
Feeds no numbered table or figure: it verifies the formulas of Proposition 4 that table_asset_pricing.py
evaluates for Table 1.  Part A is skipped if sympy is absent; the rest needs numpy and scipy.
All checks print the largest absolute error; the script stops with an AssertionError if any check fails.
Run: python3 checks_asset_pricing.py
"""
import numpy as np
from scipy.optimize import brentq

try:
    import sympy as sp
except ImportError:  # sympy is only needed for Part A
    sp = None

# ---------------------------------------------------------------- Part A
if sp is not None:
    eps, gam = sp.symbols('epsilon gamma', positive=True)
    # utility (C A^gamma)^(1-eps)/(1-eps) is homogeneous of degree (1+gamma)(1-eps) in (C, A)
    degree = (1 + gam) * (1 - eps)
    sigma_joint = 1 - degree
    print('A: 1 - degree - [eps + gamma(eps-1)] =', sp.simplify(sigma_joint - (eps + gam * (eps - 1))))

# ---------------------------------------------------------------- Part B
gh_x, gh_w = np.polynomial.hermite_e.hermegauss(80)   # nodes/weights for E f(Z), Z ~ N(0,1)
gh_w = gh_w / gh_w.sum()
E = lambda f: np.sum(gh_w * f(gh_x))


def delta_of(gamma, theta):
    """Effective curvature Delta of the isoelastic example with a linear aspiration benefit."""
    return (1 + theta - gamma) / theta


def a_closed(cbar, gamma, theta, kappa, sigma):
    """log a' from E_t U_a = -kappa when log c' ~ N(cbar, sigma^2)."""
    logEc = (1 - gamma) * cbar + 0.5 * (1 - gamma) ** 2 * sigma ** 2
    return (np.log(kappa * (gamma - 1) / (1 + theta)) - logEc) / theta


def a_numeric(cbar, gamma, theta, kappa, sigma):
    """log a' from the same first-order condition, solved numerically by bracketing the closed form."""
    def foc(loga):
        Ua = lambda z: (1 + theta) / (1 - gamma) * np.exp((1 - gamma) * (cbar + sigma * z) + theta * loga)
        return E(Ua) + kappa
    guess = a_closed(cbar, gamma, theta, kappa, sigma)
    return brentq(foc, guess - 1, guess + 1)


gamma, theta, kappa, beta, mu, sigma = 10.0, 18.0, 0.7, 0.99, 0.0183, 0.0357
Delta = delta_of(gamma, theta)
err = 0.0
for cbar in (-0.3, 0.0, 0.4):
    err = max(err, abs(a_closed(cbar, gamma, theta, kappa, sigma) - a_numeric(cbar, gamma, theta, kappa, sigma)))
print(f'B1: aspiration FOC, max |error| = {err:.2e}')
assert err < 1e-9

# B2: SDF M' = beta U_c(c', a')/U_c(c, a), U_c = c^-gamma a^(1+theta); iid log growth mu + sigma z
#     prediction: r_f = -log beta + Delta mu - (gamma - Delta) eps_t - gamma^2 sigma^2 / 2
err = 0.0
for eps_t in (-0.05, 0.0, 0.03):
    c_prev = 0.0
    c_now = c_prev + mu + eps_t
    a_now = a_closed(c_prev + mu, gamma, theta, kappa, sigma)       # set at t-1
    a_next = a_closed(c_now + mu, gamma, theta, kappa, sigma)       # set at t
    M = lambda z: beta * np.exp(-gamma * (mu + sigma * z) + (1 + theta) * (a_next - a_now))
    rf_num = -np.log(E(M))
    rf_formula = -np.log(beta) + Delta * mu - (gamma - Delta) * eps_t - 0.5 * gamma ** 2 * sigma ** 2
    err = max(err, abs(rf_num - rf_formula))
    # SDF form: log M = log beta - gamma dc' + (gamma - Delta) dc
    lhs = (1 + theta) * (a_next - a_now)
    err = max(err, abs(lhs - (gamma - Delta) * (mu + eps_t)))
print(f'B2: SDF and risk-free rate, max |error| = {err:.2e}')
assert err < 1e-9


# B3: log risk premium of the j-period consumption strip:  gamma sigma^2 + (j-1) Delta sigma^2
def strip_log_premium_mc(j, n=2_000_000, seed=1):
    """Monte Carlo log premium on a consumption strip of maturity j, held to maturity, per period j.

    The current surprise is set to zero, so the premium is measured at the conditional mean.  Compared in
    the loop below with the closed form (gamma + (j-1) Delta) sigma^2.
    """
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((n, j))
    dc = mu + sigma * z
    # log M_{t,t+j} = j log beta - gamma sum dc_{t+i} + (gamma - Delta) sum_{i=0}^{j-1} dc_{t+i}; take eps_t = 0
    logM = j * np.log(beta) - gamma * dc.sum(1) + (gamma - Delta) * (mu + dc[:, :-1].sum(1))
    logC = dc.sum(1)
    price = np.mean(np.exp(logM + logC))
    bond = np.mean(np.exp(logM))
    # log premium = log E[C_{t+j}] - log(price/bond)  (expected return on the strip over the j-period bond)
    return np.log(np.exp(j * mu + 0.5 * j * sigma ** 2)) - np.log(price / bond)


for j in (1, 2, 5):
    mc = strip_log_premium_mc(j)
    th = (gamma + (j - 1) * Delta) * sigma ** 2
    print(f'B3: strip j={j}: Monte Carlo {mc:.5f}, formula {th:.5f}')
    assert abs(mc - th) < 3e-4


# B4: consumption claim, P_t/C_t = K exp((gamma - Delta) dc_t)
def K_of(gamma, Delta, beta, mu, sigma):
    """Constant K in the price-consumption ratio of the consumption claim, P_t/C_t = K e^{(gamma-Delta) dc_t}.

    Raises an AssertionError when the geometric term q is not below one, in which case the claim has
    infinite value and the calibration is inadmissible.
    """
    num = np.exp(np.log(beta) + (1 - gamma) * mu + 0.5 * (1 - gamma) ** 2 * sigma ** 2)
    q = np.exp(np.log(beta) + (1 - Delta) * mu + 0.5 * (1 - Delta) ** 2 * sigma ** 2)
    assert q < 1, 'consumption claim has infinite value'
    return num / (1 - q)


K = K_of(gamma, Delta, beta, mu, sigma)
err = 0.0
for eps_t in (-0.05, 0.0, 0.03):
    dc_t = mu + eps_t
    pd_t = K * np.exp((gamma - Delta) * dc_t)
    rhs = E(lambda z: beta * np.exp(-gamma * (mu + sigma * z) + (gamma - Delta) * dc_t)
            * np.exp(mu + sigma * z) * (1 + K * np.exp((gamma - Delta) * (mu + sigma * z))))
    err = max(err, abs(pd_t - rhs) / pd_t)
print(f'B4: consumption-claim recursion, max relative error = {err:.2e}')
assert err < 1e-9
print('all checks passed')
