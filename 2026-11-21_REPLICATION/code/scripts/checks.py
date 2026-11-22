"""Symbolic check of the Le Chatelier result of Section 2 of 'Ex Ante Aspirations'.

With SymPy: (i) that the curvature of the indirect felicity, Delta, written in terms of the elasticities
A, B, C and D of the direct felicity, equals A (1 - U_ca^2/(U_cc U_aa)), the form stated in Proposition 1,
and that it is also the curvature of U~(c) = U(c, a*(c)) once da/dc = -U_ca/U_aa is substituted;
(ii) the Cobb-Douglas case U = c^(1-gamma) a^(1+theta)/(1-gamma) with a linear aspiration benefit kappa a,
the isoelastic example of the paper, where a*(c), U~(c) and the curvature (1 - gamma + theta)/theta are
available in closed form, together with the sign of the Hessian and the curvature at a fixed aspiration.

Feeds no numbered table or figure: it verifies the algebra of Section 2 and of the isoelastic example.
Every printed difference should be zero (the last line prints True).  Run: python3 checks.py"""
import sympy as sp
c, a, k, g, th = sp.symbols('c a kappa gamma theta', positive=True)

# general formula: Delta = A (1 - U_ca^2/(U_cc U_aa)) equals the curvature of the indirect felicity
U = sp.Function('U')
Ucc, Uca, Uaa, Uc = sp.symbols('U_cc U_ca U_aa U_c')
A = -c * Ucc / Uc
B, C_, D = a * Uca / Uc, c * Uca / (-sp.Symbol('U_a')), -a * Uaa / (-sp.Symbol('U_a'))
Delta = (A * D - B * C_) / D
print('Delta - A(1 - Uca^2/(Ucc Uaa)) =', sp.simplify(Delta - A * (1 - Uca**2 / (Ucc * Uaa))))
# indirect felicity: Ut'(c) = U_c, Ut''(c) = U_cc + U_ca da/dc with da/dc = -U_ca/U_aa
Ut_cc = Ucc - Uca**2 / Uaa
print('-c Ut_cc/U_c - Delta =', sp.simplify(-c * Ut_cc / Uc - Delta))

# Cobb-Douglas example U = c^(1-gamma) a^(1+theta)/(1-gamma), aspiration benefit kappa a
Uex = c**(1 - g) * a**(1 + th) / (1 - g)
Ua = sp.diff(Uex, a)
astar = sp.solve(sp.Eq(Ua, -k), a)[0]
Ut = sp.simplify((Uex + k * a).subs(a, astar))
elast = sp.simplify(-c * sp.diff(Ut, c, 2) / sp.diff(Ut, c))
print('a*(c) =', astar)
print('indirect felicity =', Ut)
print('curvature of indirect felicity =', sp.factor(elast), ' expected (1-gamma+theta)/theta')
print('check:', sp.simplify(elast - (1 - g + th) / th) == 0)
jac = sp.hessian(Uex, (c, a))
print('det Hessian sign factor:', sp.factor(sp.simplify(jac.det() / (c**(-2*g) * a**(2*th)))))
# risk aversion with a fixed (ex ante aspirations) vs a flexible
print('relative risk aversion, a fixed =', sp.simplify(-c * sp.diff(Uex, c, 2) / sp.diff(Uex, c)))
