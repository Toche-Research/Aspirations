# code/scripts

The scripts that produce everything the paper quotes. Run them from `code/`, not from here, with `PYTHONPATH=.`, so
that the package `aspirations/` beside them imports. Each writes into `code/out/`; none writes outside `code/`.

| What | Scripts |
|---|---|
| Table 1 and the two figures | `table_asset_pricing.py`, `fig_saving_growth.py`, `bound_curve.py` |
| The evidence of Section 7 and of the appendix | `bound_curve.py`, `real_rate_bound.py`, `spf_aspirations_test.py`, `spf_expectations_test.py`, `spf_pinned_rate.py`, `yogo.py`, `yogo_aspirations_test.py`, `checks_micro_euler.py` |
| Checks of the derivations | `checks.py`, `checks_symbolic.py`, `checks_asset_pricing.py`, `checks_risk_wealth.py`, `checks_strip_returns.py` |
| Description of the data | `data_description.py` |

`spf_aspirations_test.py` holds the regression tools (Newey-West covariance, the `Fit` class and its Wald test) that
`spf_expectations_test.py`, `spf_pinned_rate.py` and `bound_curve.py` import; `yogo.py` holds the loader and the
estimators that `yogo_aspirations_test.py` and `data_description.py` import.
