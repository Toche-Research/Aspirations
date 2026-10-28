<!-- Describes what belongs in this directory. -->

# code

The replication package: the code that produces every table and figure in the paper and the numbers
in them, and the outputs it writes. One number in the text, the bound of 2.6, is computed by hand
from the printed series and is not produced here.

| | |
|---|---|
| `aspirations/` | the model itself, imported by the scripts; writes nothing |
| `data/` | the raw sources as downloaded: FRED series, the Michigan RCONSUM series and Yogo's (2004) files. Read, never written |
| `scripts/` | the producers, run from here with `PYTHONPATH=.` set. The four `checks_*.py` verify the derivations against closed forms and limiting cases, and the three `*_test.py` run the empirical tests reported in the appendix |
| `out/` | everything the scripts write; generated, and reproducible. `out/data/` holds the derived series, `out/figures/` and `out/tables/` what the paper includes, and the `.txt` files the check and test output |
| `sync_to_paper.sh` | copies `out/` into `tex/`; the only thing that should write there |

Run the scripts from this directory, not from `scripts/`. Nothing writes outside the repository.
After a run, `sync_to_paper.sh` puts the results where the paper will find them.
