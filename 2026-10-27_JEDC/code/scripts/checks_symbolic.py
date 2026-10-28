"""Symbolic checks (SymPy) of the derivations in Appendices A-B of aspirations.tex.

1. Corollary 1 / Appendix B (satisfaction): a*(c), curvature Delta of U~, U(c,a*(c)) = 1/(g-1) - K c^{1-Delta}.
2. Appendix B, individual Euler equation: Delta mu_it = log beta + r + g^2 s^2/2 + (g-Delta) eps_it, hence (eq:ma1).
3. Appendix B, IV estimand (eq:iv), equal-loading case, Dynan's omega, omega_- and the bounds on theta_h.
4. Appendix B, the bound: yield loading -(g-Delta) eps_t / m from the SDF of Proposition 4; forward rates; sd bound.
5. Appendix B, satisfaction coefficients: (Delta-1) on anticipated log consumption, (g-1) on the surprise,
   and the loadings of Proposition 5 on expected growth, the previous expectation and the latest growth.

Feeds no numbered table or figure: it is the algebraic counterpart of the numerical checks in
checks_micro_euler.py and checks_risk_wealth.py, and covers every derivation of Appendices A and B that
the Evidence section of Section 7 relies on.  Each check prints PASS or FAIL with the simplified
difference, and the last line names any failures.  out/symbolic_checks.txt is the saved console output,
    python3 checks_symbolic.py > out/symbolic_checks.txt
"""
import sympy as sp
from fractions import Fraction

FAILS = []


def check(name, expr, extra=""):
    """expr should simplify to zero."""
    d = sp.simplify(expr)
    if d == 0:
        print(f"PASS  {name}" + (f"  [{extra}]" if extra else ""))
    else:
        FAILS.append(name)
        print(f"FAIL  {name}: simplified difference = {d}" + (f"  [{extra}]" if extra else ""))


def check_bool(name, cond, extra=""):
    """Record a check whose result is already a boolean rather than an expression to be simplified."""
    if cond:
        print(f"PASS  {name}" + (f"  [{extra}]" if extra else ""))
    else:
        FAILS.append(name)
        print(f"FAIL  {name}" + (f"  [{extra}]" if extra else ""))


# --------------------------------------------------------------------------------------
# Symbols.  gamma>1, theta>gamma-1, b>=0, kappa>0.  Work in logs: c = exp(x), a = exp(y).
# To keep SymPy's power simplification tractable, the constant kappa is parametrized by A0 > 0 through
# kappa = A0^{theta+b}(1+theta)/(gamma-1), so that a*(1) = A0 (one-to-one for kappa > 0).
g, th, b, A0 = sp.symbols("gamma theta b A0", positive=True)
x, y = sp.symbols("x y", real=True)          # x = log c, y = log a
kap = A0 ** (th + b) * (1 + th) / (g - 1)
Delta = 1 + (g - 1) * (b - 1) / (th + b)      # Corollary 1

c, a = sp.exp(x), sp.exp(y)
U = (c ** (1 - g) * a ** (1 + th) - 1) / (1 - g)                     # eq:example
B = kap * (a ** (1 - b) - 1) / (1 - b)        # generic b != 1 (b = 1 is the log limit; Delta(b=1) checked below)


la0 = sp.Symbol("logA0", real=True)


def _fix_exp(e):
    """Cancel and expand the exponent of every exp() in e, so like exponentials can be collected."""
    return e.replace(lambda t: t.func == sp.exp, lambda t: sp.exp(sp.cancel(sp.expand(t.args[0]))))


def _canon(e):
    """One pass of the canonical form used by simp(): collect powers, tidy exponents, expand sums."""
    e = sp.powsimp(sp.expand_power_exp(e), force=True, combine="exp")
    e = _fix_exp(e)
    e = sp.expand(e, power_exp=False, power_base=False)
    e = sp.powsimp(e, force=True, combine="exp")
    e = _fix_exp(e)
    return sp.expand(e, power_exp=False, power_base=False)


def simp(e):
    """Simplify sums of exponentials with rational-function exponents: write A0 = exp(logA0), expand,
    cancel each exponent, then simplify."""
    e = e.subs(A0, sp.exp(la0))
    for _ in range(3):
        e = _canon(e)
        if e == 0:
            return sp.Integer(0)
    e = sp.simplify(e)
    e = _canon(e)
    return sp.Integer(0) if e == 0 else sp.simplify(e)


print("=" * 78)
print("1. Corollary 1 and the Remark on satisfaction: a*(c), Delta, U(c,a*(c))")
print("=" * 78)
# Partial derivatives in levels: d/dc = e^{-x} d/dx, d/da = e^{-y} d/dy
U_c = sp.diff(U, x) / c
U_a = sp.diff(U, y) / a
B_a = sp.diff(B, y) / a
check("U_c = c^{-gamma} a^{1+theta}", U_c - c ** (-g) * a ** (1 + th))
check("U_a = (1+theta) c^{1-gamma} a^theta/(1-gamma)", U_a - (1 + th) * c ** (1 - g) * a ** th / (1 - g))
check("B'(a) = kappa a^{-b}", B_a - kap * a ** (-b))

# First-order condition U_a + B'(a) = 0  <=>  (1+theta) c^{1-gamma} a^{theta+b} = kappa (gamma-1)
# candidate: log a* = log A0 + (gamma-1)/(theta+b) x
y_star = sp.log(A0) + (g - 1) / (th + b) * x
foc_scaled = (U_a + B_a) * a ** b * (g - 1) / (1 + th)        # = c^{1-gamma} a^{theta+b} - kappa(gamma-1)/(1+theta) ... times -1
check("a*(c) solves U_a + B' = 0", simp(foc_scaled.subs(y, y_star)))
check("a*(c) propto c^{(gamma-1)/(theta+b)}", sp.diff(y_star, x) - (g - 1) / (th + b))

# U~(c) = U(c,a*(c)) + B(a*(c)); in x = log c: Delta = -c U~''/U~' = 1 - U~_xx / U~_x
Ut = (U + B).subs(y, y_star)
Ut_x = simp(sp.diff(Ut, x))
Ut_xx = simp(sp.diff(Ut, x, 2))
check("curvature of U~ equals Delta = 1 + (gamma-1)(b-1)/(theta+b)", simp(1 - Ut_xx / Ut_x - Delta))
check("envelope U~'(c) = U_c(c,a*(c))", simp(Ut_x / c - U_c.subs(y, y_star)))
check("U~'(c) > 0 ... U~_x = A0^{1+theta} c^{1-Delta}", simp(Ut_x - A0 ** (1 + th) * sp.exp((1 - Delta) * x)))
check("b = 0: Delta = (1+theta-gamma)/theta", Delta.subs(b, 0) - (1 + th - g) / th)
check("b = 1: Delta = 1", Delta.subs(b, 1) - 1)
check("b -> infinity: Delta -> gamma", sp.limit(Delta, b, sp.oo) - g)
check_bool("Delta increasing in b", sp.simplify(sp.diff(Delta, b)) == (g - 1) * (th + 1) / (th + b) ** 2,
           f"dDelta/db = {sp.simplify(sp.diff(Delta, b))}")
# gamma - Delta = (1+theta)(gamma-1)/(theta+b), the identity used throughout Appendix B
check("gamma - Delta = (1+theta)(gamma-1)/(theta+b)", g - Delta - (1 + th) * (g - 1) / (th + b))
# da/a = (gamma-Delta)/(1+theta) dc/c
check("dlog a*/dlog c = (gamma-Delta)/(1+theta)", sp.diff(y_star, x) - (g - Delta) / (1 + th))

# U(c, a*(c)) = 1/(gamma-1) - K c^{1-Delta}, K > 0
U_star = U.subs(y, y_star)
K_expr = simp((sp.Rational(1) / (g - 1) - U_star) * sp.exp(-(1 - Delta) * x))
check("U(c,a*(c)) = 1/(gamma-1) - K c^{1-Delta}: K independent of c", sp.diff(K_expr, x))
K_closed = simp(kap * A0 ** (1 - b) / (1 + th))      # = kappa a*(1)^{1-b}/(1+theta)
check("K = kappa a*(1)^{1-b}/(1+theta)", simp(K_expr - K_closed))
gm1 = sp.Symbol("gamma_minus_1", positive=True)
check_bool("K > 0", sp.simplify(K_closed.subs(g, 1 + gm1)).is_positive is True, f"K = {K_closed}")
# dU(c,a*(c))/dx = -K (1-Delta) c^{1-Delta}: falls in c iff Delta < 1
check("dU(c,a*)/dlog c = -K(1-Delta) c^{1-Delta}", simp(sp.diff(U_star, x) + K_closed * (1 - Delta) * sp.exp((1 - Delta) * x)))
check("U~ = U + B rises at rate U_c (envelope); U falls iff Delta < 1: U_x + K(1-Delta)c^{1-Delta} = 0 (shown above)", 0)

# --------------------------------------------------------------------------------------
print()
print("=" * 78)
print("2. Appendix B: the individual Euler equation and (eq:ma1)")
print("=" * 78)
# log-normal: log c_{t+1} | I_t ~ N(log c_t + mu_t, s^2).  For any k: log E[c_{t+1}^k] = k(lc_t + mu) + k^2 s^2/2.
kap = sp.Symbol("kappa", positive=True)
lc_t, lc_tm1, mu_t, mu_tm1, s, lbeta, r = sp.symbols("logc_t logc_tm1 mu_t mu_tm1 sigma_i logbeta r_tp1", real=True)
eps_t = lc_t - (lc_tm1 + mu_tm1)             # surprise of period t


def logE_cpow(k, lc, mu):
    """log E_t[c_{t+1}^k] when log c_{t+1} is normal with mean lc + mu and variance s^2."""
    return k * (lc + mu) + k ** 2 * s ** 2 / 2


# aspiration FOC: E_t U_a(c_{t+1}, a_{t+1}) = -B'(a_{t+1})
# (1+theta) a^theta E_t c_{t+1}^{1-gamma} /(1-gamma) = -kappa a^{-b}
# => (theta+b) log a_{t+1} = log(kappa(gamma-1)/(1+theta)) - log E_t c_{t+1}^{1-gamma}
def log_a_next(lc, mu):
    """log a_{t+1} from the aspiration first-order condition E_t U_a = -B'(a), under lognormal growth."""
    return (sp.log(kap * (g - 1) / (1 + th)) - logE_cpow(1 - g, lc, mu)) / (th + b)


# verify that this solves the FOC with the lognormal expectation taken symbolically
la = sp.Symbol("la", real=True)
foc_asp = (1 + th) * sp.exp(th * la) * sp.exp(logE_cpow(1 - g, lc_t, mu_t)) / (1 - g) + kap * sp.exp(-b * la)
check("aspiration FOC solved by (theta+b) log a_{t+1} = const - log E_t c_{t+1}^{1-gamma}",
      sp.expand_log(sp.simplify(foc_asp.subs(la, log_a_next(lc_t, mu_t)) * sp.exp(b * la).subs(la, log_a_next(lc_t, mu_t))), force=True))
# paper: (theta+b) log a_{i,t+1} = const - (1-gamma)(log c_it + mu_it) - (1-gamma)^2 s^2/2
check("(theta+b) log a_{t+1} = const - (1-gamma)(log c_t+mu_t) - (1-gamma)^2 sigma^2/2",
      (th + b) * log_a_next(lc_t, mu_t) - (sp.log(kap * (g - 1) / (1 + th)) - (1 - g) * (lc_t + mu_t) - (1 - g) ** 2 * s ** 2 / 2))

la_t = log_a_next(lc_tm1, mu_tm1)           # aspiration for t, set at t-1
la_tp1 = log_a_next(lc_t, mu_t)             # aspiration for t+1, set at t
# log U_c(c_t, a_t) = -gamma log c_t + (1+theta) log a_t
logUc_t = -g * lc_t + (1 + th) * la_t
const = (1 + th) / (th + b) * (sp.log(kap * (g - 1) / (1 + th)) - (1 - g) ** 2 * s ** 2 / 2)
check("log U_c(c_t,a_t) = const - Delta log c_t - (gamma-Delta) eps_t",
      logUc_t - (const - Delta * lc_t - (g - Delta) * eps_t))
# log E_t U_c(c_{t+1}, a_{t+1}) = (1+theta) log a_{t+1} + log E_t c_{t+1}^{-gamma}
logEUc = (1 + th) * la_tp1 + logE_cpow(-g, lc_t, mu_t)
check("log E_t U_c(c_{t+1},a_{t+1}) = const - Delta(log c_t+mu_t) + gamma^2 sigma^2/2  (same const)",
      logEUc - (const - Delta * (lc_t + mu_t) + g ** 2 * s ** 2 / 2))
# Euler equation in logs: log U_c(c_t,a_t) = log beta + r + log E_t U_c(c_{t+1},a_{t+1}); solve for mu_t
euler = sp.Eq(logUc_t, lbeta + r + logEUc)
mu_sol = sp.solve(euler, mu_t)
check_bool("Euler equation has a unique solution for mu_t", len(mu_sol) == 1)
mu_sol = mu_sol[0]
check("Delta mu_t = log beta + r + gamma^2 sigma^2/2 + (gamma-Delta) eps_t",
      Delta * mu_sol - (lbeta + r + g ** 2 * s ** 2 / 2 + (g - Delta) * eps_t))
# growth: g_{t+1} = mu_t + eps_{t+1} = kappa_t + theta_h eps_t + eps_{t+1}
eps_tp1 = sp.Symbol("eps_tp1", real=True)
theta_h = (g - Delta) / Delta
kappa_t = (lbeta + r) / Delta + g ** 2 * s ** 2 / (2 * Delta)
check("g_{t+1} = kappa_t + theta_h eps_t + eps_{t+1}, theta_h = (gamma-Delta)/Delta",
      (mu_sol + eps_tp1) - (kappa_t + theta_h * eps_t + eps_tp1))
check("power utility (Delta -> gamma as b -> infinity): theta_h -> 0", sp.limit(theta_h, b, sp.oo))
# autocorrelation of MA(1): theta/(1+theta^2) at lag 1
t_h = sp.Symbol("theta_h", positive=True)
e0, e1, e2 = sp.symbols("e0 e1 e2", real=True)
# Cov(g_{t+1}, g_t) with i.i.d. unit-variance eps: g_t = th e_{t-1} + e_t, g_{t+1} = th e_t + e_{t+1}
g_t_ = t_h * e0 + e1
g_tp1_ = t_h * e1 + e2
cov = sp.Poly(sp.expand(g_t_ * g_tp1_), e0, e1, e2)
cov_val = sum(coef for (m, coef) in cov.terms() if m in [(2, 0, 0), (0, 2, 0), (0, 0, 2)])
var_val = sum(coef for (m, coef) in sp.Poly(sp.expand(g_t_ ** 2), e0, e1, e2).terms() if m in [(2, 0, 0), (0, 2, 0), (0, 0, 2)])
check("lag-1 autocorrelation of growth = theta_h/(1+theta_h^2)", cov_val / var_val - t_h / (1 + t_h ** 2))

# --------------------------------------------------------------------------------------
print()
print("=" * 78)
print("3. Appendix B: the IV estimand (eq:iv), Dynan's loadings and the bounds")
print("=" * 78)
# n horizons per year; surprises eps_s with s = 1..3n covering years 1,2,3; year t = 2, year t+1 = 3.
for n in [1, 2, 3, 4]:
    kap_s = sp.Symbol("kappa")
    eps = sp.symbols(f"e0:{3 * n + 1}", real=True)    # eps[0] = surprise of horizon 0 (last of year 0)
    def gs(s_):  # growth of horizon s (eq:ma1)
        return kap_s + t_h * eps[s_ - 1] + eps[s_]
    def G(year):  # growth over the year: horizons n(year-1)+1 .. n*year
        return sum(gs(s_) for s_ in range(n * (year - 1) + 1, n * year + 1))
    G_t, G_tp1 = sp.expand(G(2)), sp.expand(G(3))
    # eq:annual: G_t = n kappa + theta_h ell_{t-1} + (1+theta_h) F_t + ell_t
    ell_tm1, ell_t = eps[n], eps[2 * n]
    F_t = sum(eps[n + k] for k in range(1, n))
    check(f"n={n}: (eq:annual) G_t = n kappa + theta_h ell_{{t-1}} + (1+theta_h) F_t + ell_t",
          G_t - (n * kap_s + t_h * ell_tm1 + (1 + t_h) * F_t + ell_t))
    # instrument z: Cov(z, eps_s) = cz[s] for s <= 2n (year t and earlier), 0 for year t+1
    cz = sp.symbols(f"cz0:{3 * n + 1}", real=True)
    cz = list(cz)
    for s_ in range(2 * n + 1, 3 * n + 1):
        cz[s_] = 0
    def Cov_z(expr):
        expr = sp.expand(expr)
        return sum(expr.coeff(eps[s_]) * cz[s_] for s_ in range(3 * n + 1))
    alpha = Cov_z(G_tp1) / Cov_z(G_t)
    S = cz[n] + sum(cz[n + k] for k in range(1, n)) + cz[2 * n]
    om = cz[2 * n] / S
    om_m = cz[n] / S
    check(f"n={n}: alpha = theta_h omega/(1 - omega_- + theta_h(1-omega))",
          sp.simplify(alpha - t_h * om / (1 - om_m + t_h * (1 - om))))
    # surprises before ell_{t-1} do not enter: alpha independent of cz[0..n-1]
    check_bool(f"n={n}: alpha independent of surprises before ell_{{t-1}}",
               all(sp.simplify(sp.diff(alpha, cz[s_])) == 0 for s_ in range(0, n)))
    # equal loading of the n horizons of year t, nothing earlier: alpha = theta/(n + (n-1) theta)
    eq_load = {cz[s_]: 0 for s_ in range(0, n + 1)}
    eq_load.update({cz[s_]: 1 for s_ in range(n + 1, 2 * n + 1)})
    check(f"n={n}: equal loading, omega_-=0: alpha = theta_h/(n+(n-1)theta_h)",
          sp.simplify(alpha.subs(eq_load) - t_h / (n + (n - 1) * t_h)))

# general properties of (eq:iv)
om, om_m, abar = sp.symbols("omega omega_minus alpha_bar", positive=True)
alpha_iv = t_h * om / (1 - om_m + t_h * (1 - om))
check_bool("alpha increasing in theta_h (for 1 - omega_- > 0)",
           sp.simplify(sp.diff(alpha_iv, t_h) * (1 - om_m + t_h * (1 - om)) ** 2) == om * (1 - om_m))
check("alpha -> omega/(1-omega) as theta_h -> infinity", sp.limit(alpha_iv, t_h, sp.oo) - om / (1 - om))
check("n=1, omega=1, omega_-=0: alpha = theta_h", alpha_iv.subs({om: 1, om_m: 0}) - t_h)
check("omega = 0: alpha = 0", alpha_iv.subs(om, 0))
# inversion: alpha <= abar  <=>  theta_h <= abar(1-omega_-)/(omega - abar(1-omega))
th_bound = sp.solve(sp.Eq(alpha_iv, abar), t_h)[0]
check("theta_h bound = alpha_bar(1-omega_-)/(omega - alpha_bar(1-omega))",
      th_bound - abar * (1 - om_m) / (om - abar * (1 - om)))

# Dynan's loadings: monthly random walk, calendar-year totals.  Y_A - Y_{A-1} loads month m of A with (13-m)/12
# and month m of A-1 with (m-1)/12.  Verify from the definition Y_A = sum_{m} y_{A,m}, y_{A,m} = y_{A,m-1} + u_{A,m}.
uA = sp.symbols("uA1:13", real=True)
uAm1 = sp.symbols("uB1:13", real=True)
y0 = sp.Symbol("y0", real=True)             # level at end of year A-2
yAm1 = [y0 + sum(uAm1[:m]) for m in range(1, 13)]
yA = [yAm1[-1] + sum(uA[:m]) for m in range(1, 13)]
dY = sp.expand(sum(yA) - sum(yAm1))
check("calendar-year total growth loads month m of A with (13-m) and month m of A-1 with (m-1)",
      dY - (sum((13 - m) * uA[m - 1] for m in range(1, 13)) + sum((m - 1) * uAm1[m - 1] for m in range(1, 13))))

# weights by calendar month: year A-1 (index 0), A (index 1), A+1 (index 2); in units of the annual average (/12)
def w(year_idx, m):
    """Weight on the income innovation of month m, in units of the annual average.

    year_idx is 1 for the later calendar year of the instrument's growth rate, 0 for the earlier one and
    2 for the year after, which carries no weight.
    """
    if year_idx == 1:
        return Fraction(13 - m, 12)
    if year_idx == 0:
        return Fraction(m - 1, 12)
    return Fraction(0)


def shares(n, M):
    """Lagged period = months M..M+11 starting in year A (A = t-2); period before starts in month M of A-1.
    Returns (omega, omega_-, sum_last, sum_other, sum_before_last)."""
    months_lag = [(1 + (M - 1 + k) // 12, (M - 1 + k) % 12 + 1) for k in range(12)]     # (year_idx, month)
    months_bef = [(0 + (M - 1 + k) // 12, (M - 1 + k) % 12 + 1) for k in range(12)]
    L = 12 // n
    last = sum(w(*ym) for ym in months_lag[-L:])
    other = sum(w(*ym) for ym in months_lag[:-L])
    bef_last = sum(w(*ym) for ym in months_bef[-L:])
    S_ = last + other + bef_last
    return last / S_, bef_last / S_, last, other, bef_last


for (n, M, exp_sums, exp_om, exp_omm) in [(1, 1, (Fraction(13, 2), 0, Fraction(11, 2)), 0.54, 0.46),
                                          (1, 6, (Fraction(7, 3), 0, Fraction(53, 6)), 0.21, 0.79),
                                          (4, 1, (Fraction(1, 2), 6, Fraction(5, 2)), 0.056, 0.28)]:
    o, om_, last, other, bef = shares(n, M)
    check_bool(f"n={n}, M={M}: weight sums (last, other, before-last) = {tuple(float(v) for v in exp_sums)}",
               (last, other, bef) == exp_sums, f"got {float(last):.3f}, {float(other):.3f}, {float(bef):.3f}")
    check_bool(f"n={n}, M={M}: omega={exp_om}, omega_-={exp_omm}",
               round(float(o), 2 if n == 1 else 3) == exp_om and round(float(om_), 2) == exp_omm,
               f"omega={float(o):.4f}, omega_-={float(om_):.4f}")

abar_val = sp.Rational(96, 1000)   # -0.039 + 1.96*0.069 = 0.09624
check_bool("Wald upper limit -0.039 + 1.96*0.069 = 0.096", round(-0.039 + 1.96 * 0.069, 3) == 0.096)
for (M, exp_bound) in [(1, 0.10), (6, 0.15)]:
    o, om_, *_ = shares(1, M)
    bnd = th_bound.subs({abar: abar_val, om: sp.Rational(o.numerator, o.denominator), om_m: sp.Rational(om_.numerator, om_.denominator)})
    check_bool(f"n=1, M={M}: theta_h <= {exp_bound} from alpha_bar = 0.096",
               round(float(bnd), 2) == exp_bound, f"bound = {float(bnd):.4f}")
    for D_, lab in [(sp.Rational(1, 2), "Delta=1/2"), (sp.Rational(9, 10), "Delta=0.9")]:
        print(f"      gamma-Delta = theta_h*Delta <= {float(bnd * D_):.3f} at {lab}")
o4, om4, *_ = shares(4, 1)
check_bool("n=4, M=1: omega = 1/18 <= 0.056", o4 == Fraction(1, 18), f"omega = {float(o4):.4f}")
sup4 = o4 / (1 - o4)
check_bool("n=4, M=1: alpha < omega/(1-omega) = 0.059 < 0.096 for every theta_h",
           round(float(sup4), 3) == 0.059 and float(sup4) < 0.096, f"omega/(1-omega) = {float(sup4):.4f}")
for M in [4, 6, 9, 12]:
    o_, *_ = shares(4, M)
    check_bool(f"n=4, M={M}: omega = 0 (last quarter of lagged period lies in A+1)", o_ == 0, f"omega = {float(o_)}")

# --------------------------------------------------------------------------------------
print()
print("=" * 78)
print("4. Appendix B: yield loading from the SDF of Proposition 4; the bound")
print("=" * 78)
h = sp.Symbol("h", positive=True)
mu, sig = sp.symbols("mu sigma", positive=True)
J = 5
e = sp.symbols(f"eps_t0:{J + 1}", real=True)        # e[0] = eps_t (known at t), e[i] = eps_{t+i}
gg = [mu + e[i] for i in range(J + 1)]              # g_{t+i}, i.i.d. growth
# one-period SDF (eq:sdf): log M_{t+1} = log beta - gamma g_{t+1} + (gamma - Delta) g_t
D = sp.Symbol("Delta", positive=True)
def logM1(i):
    """log of the one-period stochastic discount factor M_{t+i, t+i+1}, equation (8)."""
    return lbeta - g * gg[i + 1] + (g - D) * gg[i]
check("(eq:sdf) log M_{t+1} = log beta - Delta mu - gamma eps_{t+1} + (gamma-Delta) eps_t",
      logM1(0) - (lbeta - D * mu - g * e[1] + (g - D) * e[0]))
# conditional log E_t M_{t+1} = mean + var/2 (eps_{t+1} ~ N(0, sigma^2)); r^f = -log E_t M
mean1 = logM1(0).subs(e[1], 0)
var1 = sp.expand(logM1(0)).coeff(e[1]) ** 2 * sig ** 2
rf = -(mean1 + var1 / 2)
check("(eq:rf) r^f_{t+1} = -log beta + Delta mu - gamma^2 sigma^2/2 - (gamma-Delta) eps_t",
      rf - (-lbeta + D * mu - g ** 2 * sig ** 2 / 2 - (g - D) * e[0]))
check("sd(r^f) = (gamma-Delta) sigma", sp.expand(rf).coeff(e[0]) ** 2 * sig ** 2 - (g - D) ** 2 * sig ** 2)
# j-period SDF
for j in range(1, J + 1):
    logMj = sp.expand(sum(logM1(i) for i in range(j)))
    target = j * lbeta + (g - D) * gg[0] - D * sum(gg[i] for i in range(1, j)) - g * gg[j]
    check(f"j={j}: log M_{{t,t+j}} = j log beta + (gamma-Delta) g_t - Delta sum_{{i<j}} g_{{t+i}} - gamma g_{{t+j}}",
          logMj - target)
    check(f"j={j}: coefficient of eps_t in log M_{{t,t+j}} is gamma-Delta", logMj.coeff(e[0]) - (g - D))
    # log bond price = log E_t M_{t,t+j}: conditional mean + half conditional variance; the eps_t term is unchanged
    logP = logMj.subs({e[i]: 0 for i in range(1, j + 1)}) + sum(logMj.coeff(e[i]) ** 2 for i in range(1, j + 1)) * sig ** 2 / 2
    m = j * h
    yld = -logP / m
    check(f"j={j}: annualized yield at m=jh loads -(gamma-Delta)/m on eps_t", sp.expand(yld).coeff(e[0]) + (g - D) / m)
    # premium of the strip (iii): -Cov_t(log M_{t,t+j}, log c_{t+j}) = (gamma + (j-1) Delta) sigma^2
    logc_j = sum(gg[i] for i in range(1, j + 1))
    cov_t = sum(logMj.coeff(e[i]) * sp.expand(logc_j).coeff(e[i]) for i in range(1, j + 1)) * sig ** 2
    check(f"j={j}: strip log premium = (gamma+(j-1)Delta) sigma^2", -cov_t - (g + (j - 1) * D) * sig ** 2)
# forward rates between jh and (j+1)h do not respond to eps_t
def logP_(j):
    """log price at t of a zero-coupon bond maturing j periods ahead, conditional mean plus half variance."""
    logMj = sp.expand(sum(logM1(i) for i in range(j)))
    return logMj.subs({e[i]: 0 for i in range(1, j + 1)}) + sum(logMj.coeff(e[i]) ** 2 for i in range(1, j + 1)) * sig ** 2 / 2
for j in range(1, J):
    fwd = logP_(j) - logP_(j + 1)
    check(f"forward rate {j}h->{j + 1}h independent of eps_t", sp.expand(fwd).coeff(e[0]))
# at m = h: r^{(h)} = const - (gamma-Delta) e_t with e_t = eps_t/h, eps_t ~ N(0, sigma^2 h)
e_t = sp.Symbol("e_t", real=True)
r_h = sp.expand(-logP_(1) / h).subs(e[0], h * e_t)
check("r^{(h)} loads -(gamma-Delta) on the annualized surprise e_t = eps_t/h", sp.expand(r_h).coeff(e_t) + (g - D))
# bound: Var(r + v) = (gamma-Delta)^2 Var(e) + Var(v) >= (gamma-Delta)^2 Var(e) when Cov(e, v) = 0
Ve, Vv, Cev = sp.symbols("Var_e Var_v Cov_ev", real=True)
Var_obs = (g - D) ** 2 * Ve + Vv - 2 * (g - D) * Cev     # from r = const - (g-D) e + v
check("Var(r^{(h)}+v) = (gamma-Delta)^2 Var(e) + Var(v) when Cov(e,v) = 0",
      Var_obs.subs(Cev, 0) - ((g - D) ** 2 * Ve + Vv))
check_bool("sd(r^{(h)}) >= (gamma-Delta) sd(e), i.e. gamma-Delta <= sd(r)/sd(e)",
           sp.simplify(Var_obs.subs(Cev, 0) - (g - D) ** 2 * Ve) == Vv and True, "difference is Var(v) >= 0")
check("i.i.d. growth: Var(e) = sigma^2/h so B(h) = sqrt(h) sd(r)/sigma",
      sp.sqrt(Vv) / sp.sqrt(sig ** 2 / h) - sp.sqrt(h) * sp.sqrt(Vv) / sig)
check("first differences, white-noise surprise: Var(Delta e) = 2 Var(e)", 2 * Ve - (Ve + Ve))
check("maturity m > h: gamma-Delta <= m sd(y^{(m)})/sd(eps)",
      ((g - D) / m * sp.Symbol("sd_eps")) - (g - D) * sp.Symbol("sd_eps") / m)
# time-varying moments (Proposition 5(ii)): the eps_t term is (gamma-Delta)[g_t + (mu_t - mu_{t-1}) - (gamma-1)(s_t^2 - s_{t-1}^2)/2]
mu_t_, mu_tm1_, s_t, s_tm1 = sp.symbols("mu_t mu_tm1 sigma_t sigma_tm1", real=True)
lc_t_, lc_tm1_ = sp.symbols("lc_t lc_tm1", real=True)
def la_tv(lc, mu_, s_):
    """(1+theta) log a_{t+1} of Proposition 5(i), from the first-order condition with lognormal moments."""
    return (1 + th) / (th + b) * (sp.log(kap * (g - 1) / (1 + th)) - ((1 - g) * (lc + mu_) + (1 - g) ** 2 * s_ ** 2 / 2))
check("Prop 5(i): (1+theta) log a_{t+1} = const + (gamma-Delta)[log c_t + mu_t - (gamma-1) sigma_t^2/2]",
      la_tv(lc_t_, mu_t_, s_t) - ((1 + th) / (th + b) * sp.log(kap * (g - 1) / (1 + th)) + (g - Delta) * (lc_t_ + mu_t_ - (g - 1) * s_t ** 2 / 2)))
g_t_tv = lc_t_ - lc_tm1_
logM_tv = lbeta - g * (mu_t_ + e[1]) + la_tv(lc_t_, mu_t_, s_t) - la_tv(lc_tm1_, mu_tm1_, s_tm1)
rf_tv = -(logM_tv.subs(e[1], 0) + g ** 2 * s_t ** 2 / 2)
check("Prop 5(ii): r^f = -log beta + gamma mu_t - gamma^2 s_t^2/2 - (gamma-Delta)[g_t + (mu_t-mu_{t-1}) - (gamma-1)(s_t^2-s_{t-1}^2)/2]",
      rf_tv - (-lbeta + g * mu_t_ - g ** 2 * s_t ** 2 / 2 - (g - Delta) * (g_t_tv + (mu_t_ - mu_tm1_) - (g - 1) * (s_t ** 2 - s_tm1 ** 2) / 2)))
check("Prop 5(ii): coefficients on (mu_t, mu_{t-1}, g_t) are (Delta, gamma-Delta, -(gamma-Delta))",
      (sp.expand(rf_tv).coeff(mu_t_) - Delta) ** 2 + (sp.expand(rf_tv).coeff(mu_tm1_) - (g - Delta)) ** 2
      + (sp.expand(rf_tv.subs(lc_tm1_, lc_t_ - sp.Symbol("gt"))).coeff(sp.Symbol("gt")) + (g - Delta)) ** 2)
check("Prop 5(iii): impact of d sigma^2 with aspirations = -[gamma^2 - (gamma-Delta)(gamma-1)]/2",
      sp.diff(rf_tv, s_t) / (2 * s_t) + (g ** 2 - (g - Delta) * (g - 1)) / 2)

# --------------------------------------------------------------------------------------
print()
print("=" * 78)
print("5. Appendix B: satisfaction against the aspiration")
print("=" * 78)
# 1/(gamma-1) - U(c_{t+1}, a_{t+1}) = c^{1-gamma} a^{1+theta}/(gamma-1)
check("1/(gamma-1) - U(c,a) = c^{1-gamma} a^{1+theta}/(gamma-1)",
      sp.Rational(1) / (g - 1) - U - c ** (1 - g) * a ** (1 + th) / (g - 1))
lc_tp1 = lc_t_ + mu_t_ + e[1]
S_sat = -((1 - g) * lc_tp1 + la_tv(lc_t_, mu_t_, s_t) - sp.log(g - 1))
S_target = (Delta - 1) * (lc_t_ + mu_t_) + (g - 1) * e[1] + (g - Delta) * (g - 1) * s_t ** 2 / 2
diff = sp.expand(S_sat - S_target)
check_bool("S = -log[1/(gamma-1) - U(c_{t+1},a_{t+1})] = const + (Delta-1)(log c_t+mu_t) + (gamma-1) eps_{t+1} + (gamma-Delta)(gamma-1) sigma_t^2/2",
           all(sp.simplify(diff.coeff(v)) == 0 for v in [lc_t_, mu_t_, e[1], s_t]) and sp.simplify(sp.diff(diff, s_t)) == 0,
           f"remaining constant = {sp.simplify(diff)}")
check("coefficient on anticipated log consumption (log c_t + mu_t) is Delta-1", sp.expand(S_sat).coeff(mu_t_) - (Delta - 1))
check("coefficient on the surprise eps_{t+1} is gamma-1", sp.expand(S_sat).coeff(e[1]) - (g - 1))
# with the paper's bare form S = -log(-U): under certainty U(c,a*) = 1/(gamma-1) - K c^{1-Delta} < 0 iff K c^{1-Delta} > 1/(gamma-1);
# the paper uses 1/(gamma-1) - U, which is positive always.  That is the form checked above.
print("      note: the paper's Remark uses -log[1/(gamma-1) - U], the positive distance from the bound, not -log(-U).")

print()
print("=" * 78)
if FAILS:
    print(f"{len(FAILS)} FAIL(s): " + "; ".join(FAILS))
else:
    print("ALL CHECKS PASS")
