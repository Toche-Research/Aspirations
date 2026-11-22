# Ex Ante Aspirations, Risk Aversion, and Intertemporal Substitution

Replication package for the paper of that title by Patrick Toche. The package holds the code, the data, the outputs of
the code and a working paper version of the paper, which is the paper with two tables that describe the data added to
its appendix.

## The paper and the code

The paper studies a consumer who chooses an aspiration level for future consumption one aspiration horizon in advance,
and then faces consumption risk within that horizon with the aspiration already fixed. Because the aspiration is
chosen rather than inherited, the curvature of the consumer's effective felicity, `Delta`, is lower than the curvature
at a fixed aspiration, `gamma`, by a factor that measures the complementarity between consumption and aspirations. The
separation of the two curvatures gives an elasticity of intertemporal substitution of `1/Delta` on either side of one
together with a higher curvature toward consumption risk, a lower risk-free rate and a saving rate that rises with
growth. Its cost is a testable restriction: the risk-free rate falls after a positive consumption surprise, with
loading `-(gamma - Delta)`, which bounds the gap between the two curvatures.

The code does three things. It evaluates the model's closed forms to produce the table and the two figures of the
paper. It verifies the derivations of the appendix, symbolically with SymPy and numerically by quadrature, Monte Carlo
and simulation. It computes the empirical results of the Evidence section: a bound on `gamma - Delta` from the
volatility of ex ante real interest rates, Euler equation tests on survey expectations and on the international data
of Yogo (2004), and the individual Euler equation in the design of Dynan (2000). Every random seed is set in the
script that uses it, so every number is reproducible. One number in the text, the bound of 2.6, is computed by hand
from the printed series and is not produced by the code.

## Layout

| Path | Contents |
|---|---|
| `working_paper/aspirations_wp.pdf` | the working paper |
| `code/aspirations/` | the importable package: paths, the loader that assembles the quarterly data set, the regression tools and the style of the figures; writes nothing |
| `code/scripts/` | the scripts to run, from `code/` with `PYTHONPATH=.`; each writes into `code/out/` |
| `code/data/` | the data, one folder per source: `spf/`, `fred/`, `yogo2004/`; `README.md` (sources, citations, terms) and `DATA_DESCRIPTION.md` (every series) |
| `code/out/` | everything the scripts write: `tables/`, `figures/`, `data/` (CSV), `data_description/` and the text reports |

## Data

All the data used is in `code/data/` and ships with the package: nothing has to be downloaded and no preparation is
needed. `code/data/README.md` gives, for each source, the files, what the code reads from them, the citation, the
terms and the notice to be kept with any copy.

- **Survey of Professional Forecasters and real-time data** (`spf/`): four workbooks of the Federal Reserve Bank of
  Philadelphia. `pandas` needs `openpyxl` to read them.
- **FRED** (`fred/`): seven series. Four are produced by United States government agencies. Three are copyrighted by
  their producers and are reproduced for replication only: the two implied-volatility indices of the Chicago Board
  Options Exchange and the expected inflation of the University of Michigan's Surveys of Consumers.
- **International data of Yogo (2004)** (`yogo2004/`): fourteen country files, the data of Campbell (2003) as
  distributed by Yogo. His GAUSS programs, which `scripts/yogo.py` follows, are not shipped.

The result of Dynan (2000) that the paper discusses is taken from her published table. No micro panel is used, and no
restricted data are needed.

## Requirements

Python 3 with the libraries of `code/requirements.txt` (`pip install -r requirements.txt`). The shipped outputs were
produced with Python 3.13 and these versions:

| Library | Version | Needed by |
|---|---|---|
| numpy | 2.5.3 | every script |
| scipy | 1.18.1 | `checks_asset_pricing.py`, `checks_risk_wealth.py`, `spf_aspirations_test.py` and the scripts that import it, `yogo_aspirations_test.py` |
| pandas | 3.0.5 | the scripts that read the survey workbooks and the FRED files; `data_description.py` |
| matplotlib | 3.11.2 | `fig_saving_growth.py`, `bound_curve.py` |
| sympy | 1.14.0 | `checks.py`, `checks_symbolic.py`, part of `checks_asset_pricing.py` |
| openpyxl | 3.1.5 | read by `pandas` for the survey workbooks; not imported directly |

The two figure scripts typeset their labels with LaTeX and need `pdflatex`. `statsmodels` is not used: the Newey-West
covariance and the Wald test are written out in `numpy` in `scripts/spf_aspirations_test.py`.

## How to run

From `code/`. The whole run takes about one minute.

```
cd code
export PYTHONPATH=.

# Table 1 and the two figures
python3 scripts/table_asset_pricing.py
python3 scripts/fig_saving_growth.py
python3 scripts/bound_curve.py

# the empirical results of Section 7
python3 scripts/spf_aspirations_test.py
python3 scripts/spf_expectations_test.py
python3 scripts/spf_pinned_rate.py
python3 scripts/yogo.py                  > out/yogo_replication.txt
python3 scripts/yogo_aspirations_test.py > out/yogo_aspirations_test.txt
python3 scripts/real_rate_bound.py

# checks of the derivations of the appendix
python3 scripts/checks_symbolic.py       > out/symbolic_checks.txt
python3 scripts/checks_micro_euler.py
python3 scripts/checks.py
python3 scripts/checks_asset_pricing.py
python3 scripts/checks_risk_wealth.py
python3 scripts/checks_strip_returns.py

# the tables that describe the data (after bound_curve.py)
python3 scripts/data_description.py
```

## What a run reproduces

All sixteen scripts were run from a clean copy of this package on 22 November 2026 (Python 3.13.18, the library
versions above). Every script ended without error; the run took 59 seconds, of which `yogo_aspirations_test.py`
took 33 and `checks_micro_euler.py` 11.

- The body of Table 1 and every text report and CSV file are identical to those shipped, with one exception: the two
  files `out/data/real_rate_bound*.csv` differ from an earlier run in the sixteenth significant digit of a few numbers.
- The two figures were compared with those of the paper as images and are identical. The PDF files themselves differ
  from one run to the next in the date they carry.

| Script | What it computes | What it writes (under `code/out/`) | Where it is used |
|---|---|---|---|
| `table_asset_pricing.py` | the asset-pricing implications of the model at the calibrations of Table 1 | `tables/tbl_asset_pricing.tex` | Table 1 |
| `fig_saving_growth.py` | the saving rate against the growth rate | `figures/fig_saving_growth.pdf` | Figure 1 |
| `bound_curve.py` | the bound on `gamma - Delta` at the quarterly and annual horizons, with its variants | `data/bound_curve.csv`, `data/bound_curve_data_quarterly.csv`, `data/bound_curve_data_annual.csv`, `figures/fig_bound_curve.pdf` | Figure 2 and the paragraph "The bound" of Section 7 |
| `real_rate_bound.py` | the annual bound from FRED data only, a check on `bound_curve.py` | `data/real_rate_bound.csv`, `data/real_rate_bound_data.csv` | Section 7 |
| `spf_aspirations_test.py` | the regressions of the real rate on survey expectations | `data/spf_test_results.csv`, `data/spf_test_data.csv` | the appendix on the evidence |
| `spf_expectations_test.py` | the regressions solved for expected growth | `data/spf_expectations_results.csv` | the appendix on the evidence |
| `spf_pinned_rate.py` | the regressions when the bill rate is above 0.25 per cent, and with a measure of risk | `data/spf_pinned_rate_results.csv` | the appendix on the evidence |
| `yogo.py` | the replication of the estimates of Yogo (2004) | `yogo_replication.txt` | the appendix on the evidence |
| `yogo_aspirations_test.py` | the tests of the model on the data of Yogo (2004) | `yogo_aspirations_test.txt` | the appendix on the evidence |
| `checks_micro_euler.py` | the individual Euler equation in the design of Dynan (2000), by simulation | `micro_euler_dynan.txt` | Section 7 and the appendix |
| `checks.py`, `checks_symbolic.py`, `checks_asset_pricing.py`, `checks_risk_wealth.py`, `checks_strip_returns.py` | checks of the derivations | `symbolic_checks.txt`; the others print to the console | the appendix on proofs |
| `data_description.py` | for every series: source, period, observations, mean and standard deviation | `data_description/` and `data/DATA_DESCRIPTION.md` | the two data tables of the working paper |

## Licences

Code: MIT (`code/LICENSE`). Outputs: CC BY 4.0. Data: the terms of each source. See `code/LICENSES.md` and
`code/data/README.md`.
