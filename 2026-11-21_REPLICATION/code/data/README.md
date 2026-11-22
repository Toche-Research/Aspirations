# code/data

The data the paper uses. All of it is shipped here; nothing has to be downloaded.

None of these data are the author's. Each belongs to the source named, under the terms stated; the notices in this file
are to be kept with any copy. The series, tables and figures that the scripts compute from them are in `../out/`.

| Folder | Source | Shipped | Used by |
|---|---|---|---|
| `spf/` | Survey of Professional Forecasters and real-time data on consumption, Federal Reserve Bank of Philadelphia | yes | `aspirations/data.py` and the scripts that import it; `scripts/bound_curve.py` |
| `fred/` | United States interest rates, expected inflation, consumption, prices and implied volatility, retrieved from FRED | yes | `aspirations/data.py`, `scripts/bound_curve.py`, `scripts/real_rate_bound.py`, `scripts/spf_pinned_rate.py` |
| `yogo2004/` | International data of Campbell (2003), as distributed by Yogo (2004) | yes | `scripts/yogo.py`, `scripts/yogo_aspirations_test.py` |

`DATA_DESCRIPTION.md` describes every series used in the paper: source, period, observations, mean and standard
deviation. It is written by `scripts/data_description.py`.

## spf/ (shipped)

Four Excel workbooks downloaded on 17 October 2026 from the Federal Reserve Bank of Philadelphia
(https://www.philadelphiafed.org/surveys-and-data), unchanged.

| File | Contents | Read |
|---|---|---|
| `Mean_RCONSUM_Growth.xlsx` | Survey of Professional Forecasters: mean forecasts of the growth of real personal consumption expenditures, by survey quarter (`YEAR`, `QUARTER`) and horizon (`DRCONSUM2` the current quarter, `DRCONSUM3` to `DRCONSUM6` the next four) | `DRCONSUM2` to `DRCONSUM6` |
| `Mean_PGDP_Growth.xlsx` | The same for the GDP price index (`DPGDP2`, `DPGDP3`, ...) | `DPGDP2`, `DPGDP3` |
| `Dispersion_RCONSUM.xlsx` | Cross-sectional dispersion of the forecasts of real consumption growth | sheet `D2`, the current quarter |
| `rcon_first_second_third.xlsx` | Real-time data: first, second and third releases of the growth of real personal consumption expenditures, by quarter | sheet `DATA`, columns `First` and `Second` |

Growth rates are annualized, in per cent.

Notice. The data are those of the Federal Reserve Bank of Philadelphia and are reproduced here unchanged, for the
purpose of replication. The workbooks state no terms of use; anyone who wishes to reuse them should take them from the
Bank's site and cite: Federal Reserve Bank of Philadelphia, Survey of Professional Forecasters, and Real-Time Data Set
for Macroeconomists (first, second and third releases), retrieved 17 October 2026.

## fred/ (shipped)

Seven series of the United States, one CSV file each, downloaded from FRED, Federal Reserve Bank of St. Louis
(https://fred.stlouisfed.org), on 17 October 2026. The file name is the FRED series identifier. Each file has a
column `observation_date` and a column named after the series.

| File | Series | Original source | Rights, as tagged by FRED |
|---|---|---|---|
| `TB3MS.csv` | 3-month Treasury bill, secondary market rate, monthly | Board of Governors of the Federal Reserve System | public domain |
| `GS1.csv` | Market yield on U.S. Treasury securities at 1-year constant maturity, monthly | Board of Governors of the Federal Reserve System | public domain |
| `A794RX0Q048SBEA.csv` | Real personal consumption expenditures per capita, quarterly | U.S. Bureau of Economic Analysis | public domain |
| `A191RI1Q225SBEA.csv` | Gross domestic product: implicit price deflator, percentage change from the preceding period, annual rate, quarterly | U.S. Bureau of Economic Analysis | public domain |
| `MICH.csv` | University of Michigan: inflation expectation, median expected price change over the next 12 months, monthly | University of Michigan, Surveys of Consumers | copyrighted |
| `VIXCLS.csv` | CBOE Volatility Index: VIX, daily close | Chicago Board Options Exchange | copyrighted: citation required |
| `VXOCLS.csv` | CBOE S&P 100 Volatility Index: VXO, daily close | Chicago Board Options Exchange | copyrighted |

Notices.

- The first four series are produced by agencies of the United States government and are in the public domain. They
  were accessed via FRED, whose terms of use ask users to cite the source and to note that the data were accessed via
  FRED (https://fred.stlouisfed.org/legal/). Citation: Board of Governors of the Federal Reserve System and U.S. Bureau
  of Economic Analysis, series named above, retrieved from FRED, Federal Reserve Bank of St. Louis, 17 October 2026.
- `VIXCLS.csv` and `VXOCLS.csv`: Copyright, Chicago Board Options Exchange, Inc. Reprinted with permission (the notice
  FRED carries for these series). Citation: Chicago Board Options Exchange, CBOE Volatility Index: VIX [VIXCLS] and
  CBOE S&P 100 Volatility Index: VXO [VXOCLS], retrieved from FRED, Federal Reserve Bank of St. Louis;
  https://fred.stlouisfed.org/series/VIXCLS and https://fred.stlouisfed.org/series/VXOCLS, 17 October 2026.
- `MICH.csv`: © University of Michigan, Surveys of Consumers. Citation: University of Michigan, University of Michigan:
  Inflation Expectation [MICH], retrieved from FRED, Federal Reserve Bank of St. Louis;
  https://fred.stlouisfed.org/series/MICH, 17 October 2026.
- The three copyrighted series are reproduced here unchanged, for the purpose of replication only. No licence is
  granted here over them; anyone who wishes to reuse them should take them from FRED and follow the notes of each
  series there.
- The series titles in the table are those FRED gives; check them against FRED before relying on them.

## yogo2004/ (shipped)

Fourteen plain-text files, one per country and frequency: `<CTRY>Q.txt` quarterly (Australia `AUL`, Canada `CAN`,
France `FR`, Germany `GER`, Italy `ITA`, Japan `JAP`, the Netherlands `NTH`, Sweden `SWD`, Switzerland `SWT`, the
United Kingdom `UK`, the United States `USA`) and `<CTRY>Y.txt` annual (Sweden, the United Kingdom, the United States).
They are the data of Campbell (2003) as distributed by Motohiro Yogo with Yogo (2004), in the archive `EIS_Data.zip`
(originally at http://www.nber.org/~myogo/data/EIS_Data.zip), unpacked and unchanged; the files are dated 11 November
2002.

Columns: `DATE` (year, and quarter after the point), `r` log stock return, `dp` log dividend-price ratio, `rf` log
nominal interest rate, `inf` log inflation, `dc` log consumption growth, `rr` log real stock return, `rrf` log real
interest rate, `z1` to `z4` the instruments (the log dividend-price ratio, the nominal interest rate, inflation and consumption
growth, lagged twice). A missing value is written `.`.

Notice. The data are the work of John Y. Campbell and Motohiro Yogo, built from sources that those papers name, and are
reproduced here unchanged, for the purpose of replication only. The archive states no terms of use; no licence is
granted here over them, and anyone who wishes to reuse them should obtain them from the authors and cite: Yogo,
Motohiro (2004), "Estimating the Elasticity of Intertemporal Substitution When Instruments Are Weak", Review of
Economics and Statistics, 86(3), 797–810; and Campbell, John Y. (2003), "Consumption-Based Asset Pricing", in
Handbook of the Economics of Finance, edited by George M. Constantinides, Milton Harris and René M. Stulz, Amsterdam:
Elsevier (the source Yogo names for the data).
