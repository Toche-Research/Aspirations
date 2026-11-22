"""Strip premia in the aspirations economy (Proposition 4): hold-to-maturity versus one-period holding returns.

SDF: log M' = log beta - Delta*mu - gamma*e' + (gamma - Delta)*e, consumption growth g' = mu + e', e' ~ N(0, sigma^2) i.i.d.
Closed forms (log premia over the bond of the same holding period):
  hold to maturity, per year, maturity j:  (gamma + (j-1)*Delta) * sigma^2 / j   -> falls from gamma*sigma^2 to Delta*sigma^2
  one-period holding return, maturity j:   gamma*sigma^2 for j = 1,  gamma*(1 + gamma - Delta)*sigma^2 for every j >= 2
(the price-consumption ratio of a strip loads (gamma - Delta) on the current shock, so a strip that does not pay next
period gains when consumption surprises upward). Checked here by Monte Carlo, with the pricing condition E[M R] = 1.
Calibration of the paper's table: mu = 1.83%, sigma = 3.57% (Mehra-Prescott 1889-1978), beta = 0.97; gamma = 2, Delta = 0.6.

Feeds no numbered table or figure directly: it verifies part (iii) of Proposition 4, the term structure of
strip premia that lies behind the two zero-coupon claim columns of Table 1 (computed in
table_asset_pricing.py).  Everything is printed; nothing is written to disk.  Run: python3 checks_strip_returns.py"""
import numpy as np
b, mu, s, g, D = 0.97, 0.0183, 0.0357, 2.0, 0.6
N = 2_000_000
rng = np.random.default_rng(0)

def logPc(j, e):
    """log(P^j_t / c_t) given the current shock e_t (closed form)."""
    return j*(np.log(b) - D*mu + mu) + (g - D)*e + 0.5*((1 - g)**2 + (j - 1)*(1 - D)**2)*s**2

def one_period_log_return(j, e0, e1):
    """Log return over one period on a strip of maturity j, given shocks e_t = e0 and e_{t+1} = e1.

    A one-period strip pays consumption next period; a longer strip is sold on at maturity j - 1.
    """
    return (mu + e1 - logPc(1, e0)) if j == 1 else (logPc(j - 1, e1) + mu + e1 - logPc(j, e0))

e0, e1 = rng.normal(0, s, N), rng.normal(0, s, N)
logM = np.log(b) - D*mu - g*e1 + (g - D)*e0
print("pricing check E[M R] (should be 1):", [round(float(np.mean(np.exp(logM + one_period_log_return(j, e0, e1)))), 4) for j in (1, 2, 5, 10)])
e0 = np.zeros(N)                                 # conditional on e_t = 0
lRf = -(np.log(b) - D*mu + 0.5*g**2*s**2)
print("one-period log premium, percent: maturity, Monte Carlo, closed form")
for j in (1, 2, 5, 10):
    lR = one_period_log_return(j, e0, e1)
    print(j, round(100*(lR.mean() - lRf + 0.5*lR.var()), 3), round(100*(g*s*s if j == 1 else g*(1 + g - D)*s*s), 3))
print("hold-to-maturity log premium per year, percent:", {j: round(100*(g + (j - 1)*D)*s*s/j, 3) for j in (1, 2, 5, 10, 50)})
