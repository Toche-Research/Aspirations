# code/out

Everything the scripts write. The whole tree is generated and is reproduced by running the scripts again; nothing in
it should be edited by hand. The outputs are shipped so that a reader can compare a run with them.

| Path | What it holds | Written by |
|---|---|---|
| `tables/tbl_asset_pricing.tex` | the body of Table 1 | `table_asset_pricing.py` |
| `figures/fig_saving_growth.pdf`, `figures/fig_bound_curve.pdf` | Figures 1 and 2 | `fig_saving_growth.py`, `bound_curve.py` |
| `data/bound_curve*.csv` | the bounds, one row per variant, and the assembled quarterly and annual series | `bound_curve.py` |
| `data/real_rate_bound*.csv` | the annual bound from FRED data only, and its series | `real_rate_bound.py` |
| `data/spf_test_*.csv`, `data/spf_expectations_results.csv`, `data/spf_pinned_rate_results.csv` | the regressions on survey expectations, one row per specification, and the assembled quarterly series | `spf_aspirations_test.py`, `spf_expectations_test.py`, `spf_pinned_rate.py` |
| `yogo_replication.txt`, `yogo_aspirations_test.txt` | the replication of Yogo (2004) and the tests on his data | `yogo.py`, `yogo_aspirations_test.py` |
| `micro_euler_dynan.txt` | the individual Euler equation in the design of Dynan (2000) | `checks_micro_euler.py` |
| `symbolic_checks.txt` | the symbolic checks of the derivations | `checks_symbolic.py` |
| `data_description/` | the tables that describe the data, the same in Markdown, and every number unrounded | `data_description.py` |

The other check scripts print to the console and write no file.
