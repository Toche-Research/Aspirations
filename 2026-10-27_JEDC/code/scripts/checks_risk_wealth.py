"""Numerical checks of the curvature, certainty-equivalent and wealth results of 'Ex Ante Aspirations'.

C1  Generalized benefit B(a) = kappa a^(1-b)/(1-b) (log for b = 1): the effective curvature
    Delta = 1 + (gamma - 1)(b - 1)/(theta + b), by numerical maximization over a and finite differences of U~.
C2  Certainty-equivalent representation: max_a E[U(c, a)] + B(a) = K * CE(c)^(1-Delta)/(1-Delta) with
    CE(c) = (E c^(1-gamma))^(1/(1-gamma)), for lognormal c, by quadrature.
C3  Relative risk aversion over wealth when the current aspiration is fixed and future aspirations adjust:
    RRA_W = 1/(s/gamma + (1-s)/Delta), s = c/W, against a direct computation of J(W) for a consumer who consumes c
    now and saves W - c at gross return R into a future value F(S) = K S^(1-Delta)/(1-Delta).
C4  Time-varying conditional moments (Proposition 5): the aspiration rule and the risk-free rate, by quadrature.

Feeds no numbered table or figure: C1 to C3 verify Corollary 1 and Proposition 3, which Section 3 and
Section 4 state, and C4 verifies Proposition 5, which the Euler-equation tests of Section 7 are built on.
Every check prints the formula beside the numerical value; the script stops with an AssertionError if the
largest relative error exceeds its tolerance.  Run: python3 checks_risk_wealth.py
"""
import numpy as np
from scipy.optimize import minimize_scalar, brentq

gh_x, gh_w = np.polynomial.hermite_e.hermegauss(80)
gh_w = gh_w / gh_w.sum()


def U(c, a, g, th):
    """Felicity of the isoelastic example, c^(1-gamma) a^(1+theta)/(1-gamma)."""
    return c ** (1 - g) * a ** (1 + th) / (1 - g)


def B(a, k, b):
    """Generalized aspiration benefit kappa a^(1-b)/(1-b), with the logarithmic case at b = 1."""
    return k * np.log(a) if b == 1 else k * a ** (1 - b) / (1 - b)


def a_opt(EU_coef, g, th, k, b):
    """argmax_a a^(1+th) * EU_coef/(1-g) + B(a): closed form a^(th+b) = -k (1-g) / ((1+th) EU_coef)"""
    return (-k * (1 - g) / ((1 + th) * EU_coef)) ** (1 / (th + b))


def Utilde(c, g, th, k, b):
    """Indirect felicity max_a U(c, a) + B(a) under certainty, at the optimal aspiration."""
    a = a_opt(c ** (1 - g), g, th, k, b)
    return U(c, a, g, th) + B(a, k, b)


def delta_formula(g, th, b):
    """Curvature of the indirect felicity, Delta = 1 + (gamma-1)(b-1)/(theta+b) (Corollary 1)."""
    return 1 + (g - 1) * (b - 1) / (th + b)


worst = 0.0
# C1
for g, th, b in [(2, 5, 0), (3, 4, 0.5), (5, 10, 1), (3, 4, 2), (8, 20, 5)]:
    k = 0.7
    # numerical check of the inner maximization
    c = 1.3
    res = minimize_scalar(lambda la: -(U(c, np.exp(la), g, th) + B(np.exp(la), k, b)), bounds=(-20, 20), method='bounded',
                          options={'xatol': 1e-12})
    err_a = abs(np.exp(res.x) - a_opt(c ** (1 - g), g, th, k, b)) / a_opt(c ** (1 - g), g, th, k, b)
    h = 1e-4
    up = (Utilde(c + h, g, th, k, b) - Utilde(c - h, g, th, k, b)) / (2 * h)
    upp = (Utilde(c + h, g, th, k, b) - 2 * Utilde(c, g, th, k, b) + Utilde(c - h, g, th, k, b)) / h ** 2
    D_num = -c * upp / up
    err = abs(D_num - delta_formula(g, th, b))
    worst = max(worst, err, err_a)
    print(f'C1 gamma={g} theta={th} b={b}: Delta formula {delta_formula(g, th, b):.6f}, numerical {D_num:.6f}, '
          f'argmax error {err_a:.1e}')

# C2
for g, th, b in [(2, 5, 0), (3, 4, 2), (5, 10, 1)]:
    k, mu, s = 0.7, 0.02, 0.2
    c = np.exp(mu + s * gh_x)
    EUc = np.sum(gh_w * c ** (1 - g))
    a = a_opt(EUc, g, th, k, b)
    val = a ** (1 + th) * EUc / (1 - g) + B(a, k, b)
    CE = EUc ** (1 / (1 - g))
    D = delta_formula(g, th, b)
    val_det = Utilde(CE, g, th, k, b)     # the same aspiration problem with consumption CE for sure
    err = abs(val - val_det) / abs(val)
    worst = max(worst, err)
    print(f'C2 gamma={g} theta={th} b={b}: max_a E[U]+B = {val:.8f}, U~(CE) = {val_det:.8f} (U~ is isoelastic, '
          f'curvature {D:.4f})')

# C3
for g, th, b, s_share in [(2, 5, 0, 0.03), (5, 10, 0, 0.03), (10, 18, 0, 0.03), (3, 4, 2, 0.2)]:
    k = 1.0
    D = delta_formula(g, th, b)
    R = 1.03
    W0 = 1.0
    c0 = s_share * W0
    a0 = a_opt(c0 ** (1 - g), g, th, k, b)               # aspiration set for the expected consumption
    Uc = lambda c: c ** (-g) * a0 ** (1 + th)
    # choose K so that c0 is optimal at W0: Uc(c0) = R F'(R (W0 - c0)) with F(x) = K x^(1-D)/(1-D)
    S0 = W0 - c0
    K = Uc(c0) / (R * (R * S0) ** (-D))
    F = lambda x: K * x ** (1 - D) / (1 - D)

    def J(W):
        f = lambda c: Uc(c) - R * K * (R * (W - c)) ** (-D)
        c = brentq(f, 1e-9 * W, W * (1 - 1e-12))
        return U(c, a0, g, th) + B(a0, k, b) + F(R * (W - c))

    h = 1e-3
    J1 = (J(W0 + h) - J(W0 - h)) / (2 * h)
    J2 = (J(W0 + h) - 2 * J(W0) + J(W0 - h)) / h ** 2
    rra = -W0 * J2 / J1
    formula = 1 / (s_share / g + (1 - s_share) / D)
    err = abs(rra - formula)
    worst = max(worst, err / formula)
    print(f'C3 gamma={g} Delta={D:.3f} s={s_share}: RRA over wealth {rra:.5f}, formula {formula:.5f}')

print('largest relative error', worst)
assert worst < 1e-4

# C4: time-varying conditional moments. log c_{t+1} - log c_t ~ N(mu_t, s_t^2) given date-t information.
# Claim: (1+theta) log a_{t+1} = const + (gamma - Delta)[log c_t + mu_t - (gamma-1) s_t^2/2]
#        r^f_{t+1} = -log beta + gamma mu_t - gamma^2 s_t^2/2 - (gamma - Delta)[g_t + (mu_t - mu_{t-1}) - (gamma-1)(s_t^2 - s_{t-1}^2)/2]
def check_c4():
    """Check Proposition 5 by quadrature at random date-t moments, and return the largest error.

    For each draw the aspiration is solved from its first-order condition at t-1 and at t, the risk-free
    rate is computed from the stochastic discount factor by Gauss-Hermite quadrature, and the result is
    compared with the closed form, in which the latest growth, the change in expected growth and the
    change in conditional variance all enter with coefficient (gamma - Delta).
    """
    worst = 0.0
    beta_ = 0.97
    for g, th, b in [(3.0, 4.0, 0.0), (5.0, 8.0, 2.0), (2.0, 3.0, 1.0)]:
        k = 0.8
        D = delta_formula(g, th, b)
        rng = np.random.default_rng(0)
        for _ in range(5):
            lc_prev, mu_prev, s_prev = 0.3 * rng.standard_normal(), 0.02 + 0.01 * rng.standard_normal(), 0.03 + 0.02 * rng.random()
            g_t = mu_prev + s_prev * rng.standard_normal()
            lc_t = lc_prev + g_t
            mu_t, s_t = 0.02 + 0.01 * rng.standard_normal(), 0.03 + 0.02 * rng.random()
            # aspirations from the first-order condition E U_a + B' = 0, by quadrature
            def a_of(lc, mu, s):
                EU = np.sum(gh_w * np.exp((1 - g) * (lc + mu + s * gh_x)))
                return a_opt(EU, g, th, k, b)
            a_t = a_of(lc_prev, mu_prev, s_prev)          # set at t-1 for date t
            a_n = a_of(lc_t, mu_t, s_t)                   # set at t for date t+1
            # SDF and risk-free rate by quadrature
            Uc_t = np.exp(-g * lc_t) * a_t ** (1 + th)
            M = beta_ * np.exp(-g * (lc_t + mu_t + s_t * gh_x)) * a_n ** (1 + th) / Uc_t
            rf_num = -np.log(np.sum(gh_w * M))
            rf_formula = (-np.log(beta_) + g * mu_t - 0.5 * g ** 2 * s_t ** 2
                          - (g - D) * (g_t + (mu_t - mu_prev) - 0.5 * (g - 1) * (s_t ** 2 - s_prev ** 2)))
            worst = max(worst, abs(rf_num - rf_formula))
    print(f'C4 time-varying moments: largest |r^f quadrature - formula| = {worst:.2e}')
    return worst

w4 = check_c4()
assert w4 < 1e-10
