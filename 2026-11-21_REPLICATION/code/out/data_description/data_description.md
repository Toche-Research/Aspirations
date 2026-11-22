# Data description

Written by `code/scripts/data_description.py`; do not edit. The same tables are in the appendix of the working paper.

## United States: survey expectations, interest rates and consumption

All series are in per cent a year. SPF: Survey of Professional Forecasters, Federal Reserve Bank of Philadelphia (`spf/`).
Real-time: first releases of the growth of real consumption, Federal Reserve Bank of Philadelphia (`spf/rcon_first_second_third.xlsx`).
FRED: `fred/`. Each row covers the observations of the series itself; the regressions use the quarters in which
all their variables are present. The column "Variable" is the name of the series in `code/out/data/bound_curve_data_quarterly.csv`
and `bound_curve_data_annual.csv`.

| Series | Variable | Source | Frequency | First | Last | Observations | Mean | S.d. |
|---|---|---|---|---|---|---:|---:|---:|
| Expected growth, current quarter | `mu` | SPF | quarterly | 1981Q3 | 2026Q3 | 181 | 2.40 | 3.52 |
| Expected growth, next quarter | `mu_next` | SPF | quarterly | 1981Q3 | 2026Q3 | 181 | 2.54 | 1.28 |
| Expected inflation over the bill's life | `pi_e` | SPF | quarterly | 1968Q4 | 2026Q3 | 232 | 3.39 | 2.03 |
| Dispersion of the growth forecasts | `disp` | SPF | quarterly | 1981Q3 | 2026Q3 | 181 | 1.36 | 1.68 |
| Growth of real consumption, first release | `g` | Real-time | quarterly | 1965Q4 | 2026Q3 | 242 | 2.95 | 4.47 |
| Growth surprise | `e` | SPF, Real-time | quarterly | 1981Q4 | 2026Q3 | 178 | 0.44 | 2.24 |
| Three-month Treasury bill rate | `tbill` | FRED | quarterly | 1934Q1 | 2026Q3 | 371 | 3.42 | 3.09 |
| Ex ante real bill rate | `r` | FRED, SPF | quarterly | 1968Q4 | 2026Q3 | 232 | 1.10 | 2.20 |
| Ex post real bill rate | `r_expost` | FRED | quarterly | 1947Q1 | 2026Q1 | 317 | 0.73 | 2.90 |
| Growth of real consumption per capita | `g_pce` | FRED | quarterly | 1947Q3 | 2026Q3 | 317 | 2.08 | 4.27 |
| Ex ante real one-year rate, February | `r_feb` | FRED | annual | 1978 | 2026 | 49 | 1.14 | 3.08 |
| Ex post real one-year rate, February | `r_feb_expost` | FRED | annual | 1954 | 2026 | 73 | 1.32 | 2.73 |
| Expected growth over the year, first-quarter survey | `mu_q1` | SPF | annual | 1982 | 2026 | 45 | 2.55 | 0.85 |
| Growth of real consumption over the year, first release | `g` | Real-time | annual | 1966 | 2024 | 59 | 2.93 | 1.65 |
| Growth surprise over the year | `e_q1` | SPF, Real-time | annual | 1982 | 2024 | 43 | 0.33 | 1.17 |
| Growth of real consumption per capita | `g_pce` | FRED | annual | 1948 | 2025 | 78 | 2.08 | 1.87 |

## International data of Yogo (2004)

The rows that the estimates use: observations with a missing value of any variable or instrument are dropped, as in
Yogo's programs. Log growth and log returns, in per cent a year (quarterly values are multiplied by four). The second
United States quarterly series starts in 1970, as in Yogo's short sample.

| File | Country | Frequency | First | Last | Observations | Consumption growth: mean | s.d. | Real interest rate: mean | s.d. | Real stock return: mean | s.d. |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| `USAQ.txt` | United States | quarterly | 1947Q3 | 1998Q4 | 206 | 2.0 | 2.1 | 0.9 | 3.5 | 8.1 | 31.4 |
| `AULQ.txt` | Australia | quarterly | 1970Q3 | 1998Q4 | 114 | 2.1 | 4.1 | 2.1 | 5.1 | 3.9 | 45.8 |
| `CANQ.txt` | Canada | quarterly | 1970Q3 | 1999Q1 | 115 | 2.2 | 3.9 | 2.7 | 3.8 | 5.8 | 33.7 |
| `FRQ.txt` | France | quarterly | 1970Q3 | 1998Q3 | 113 | 1.6 | 5.8 | 2.5 | 3.7 | 8.4 | 45.3 |
| `GERQ.txt` | Germany | quarterly | 1979Q1 | 1998Q3 | 79 | 1.7 | 4.8 | 3.2 | 2.3 | 10.1 | 42.3 |
| `ITAQ.txt` | Italy | quarterly | 1971Q4 | 1998Q1 | 106 | 2.3 | 3.3 | 2.4 | 5.8 | 4.4 | 54.3 |
| `JAPQ.txt` | Japan | quarterly | 1970Q3 | 1998Q4 | 114 | 3.2 | 5.1 | 1.3 | 4.6 | 4.7 | 43.4 |
| `NTHQ.txt` | Netherlands | quarterly | 1977Q3 | 1998Q4 | 86 | 1.8 | 5.0 | 3.4 | 3.1 | 14.1 | 34.7 |
| `SWDQ.txt` | Sweden | quarterly | 1970Q3 | 1999Q2 | 116 | 1.0 | 3.7 | 2.0 | 5.7 | 11.7 | 47.7 |
| `SWTQ.txt` | Switzerland | quarterly | 1976Q2 | 1998Q4 | 91 | 0.8 | 4.9 | 0.6 | 3.4 | 9.9 | 39.6 |
| `UKQ.txt` | United Kingdom | quarterly | 1970Q3 | 1999Q1 | 115 | 2.2 | 5.1 | 1.3 | 6.0 | 9.0 | 42.5 |
| `USAQ.txt` | United States (from 1970) | quarterly | 1970Q3 | 1998Q4 | 114 | 1.8 | 1.8 | 1.5 | 3.4 | 8.0 | 33.9 |
| `SWDY.txt` | Sweden | annual | 1921 | 1994 | 74 | 1.8 | 2.8 | 2.0 | 5.9 | 6.6 | 18.5 |
| `UKY.txt` | United Kingdom | annual | 1921 | 1994 | 74 | 1.5 | 2.6 | 1.4 | 5.4 | 7.7 | 22.1 |
| `USAY.txt` | United States | annual | 1891 | 1995 | 105 | 1.7 | 3.2 | 2.0 | 8.9 | 6.7 | 18.6 |
