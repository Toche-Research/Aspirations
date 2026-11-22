"""Figure 1 of the paper (Section 5, 'Saving and growth'): steady-state saving rate against growth.

Plots s = alpha (g + delta)/(delta + beta + Delta g), the saving rate on the balanced growth path of the
neoclassical economy of Proposition 6, against the growth rate g, for gamma = 3 without aspirations
(Delta = gamma) and with aspirations of strength theta = 20, 4 and 2.5 in the isoelastic example with a
linear benefit, so that Delta = (1 + theta - gamma)/theta.

Writes out/figures/fig_saving_growth.pdf and prints Delta, the saving rate at g = 0 and at g = 4 per cent,
and the slope between them, for each case.  The output path is relative, so run the script from code/."""
import os
import numpy as np
from aspirations import style  # selects the pgf backend; must precede pyplot
style.setup()
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
beta, delta, gamma, alpha = 0.05, 0.05, 3.0, 1/3
g = np.linspace(0, 0.04, 81)
fig, ax = plt.subplots(figsize=(4.6, 3.2))
for theta, ls in ((None, '-'), (20.0, '--'), (4.0, '-.'), (2.5, ':')):
    D = gamma if theta is None else (1 + theta - gamma) / theta
    lab = rf'no aspirations ($\Delta=\gamma={gamma:g}$)' if theta is None else rf'$\theta={theta:g}$ ($\Delta={D:.2f}$)'
    ax.plot(100 * g, 100 * alpha * (g + delta) / (delta + beta + D * g), ls, color='k', lw=1.2, label=lab)
ax.set_xlabel('growth rate $g$ (%)'); ax.set_ylabel('saving rate $s$ (%)')
ax.legend(frameon=False, fontsize=style.TICK)
fig.tight_layout(); FIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out', 'figures')
os.makedirs(FIG, exist_ok=True)
fig.savefig(os.path.join(FIG, 'fig_saving_growth.pdf'), bbox_inches="tight")
for theta in (None, 20.0, 4.0, 2.5):
    D = gamma if theta is None else (1 + theta - gamma) / theta
    s0, s4 = alpha * delta / (delta + beta), alpha * (0.04 + delta) / (delta + beta + D * 0.04)
    print(f'theta={theta}: Delta={D:.3f}, s(0)={s0:.4f}, s(4%)={s4:.4f}, slope={(s4 - s0) / 0.04:.3f}')
