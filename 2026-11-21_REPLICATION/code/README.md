# code

The code, the data and the outputs of the replication package. The README one level up says what the package is,
what it needs, how to run it and what a run reproduces.

| | |
|---|---|
| `aspirations/` | the importable package: paths, the loader that assembles the quarterly data set, the regression tools, the style of the figures; writes nothing |
| `scripts/` | the Python scripts, run from this folder with `PYTHONPATH=.`; see `scripts/README.md` |
| `data/` | the data, one folder per source; `data/README.md` gives sources, citations and terms, and `data/DATA_DESCRIPTION.md` describes every series |
| `out/` | everything the scripts write; see `out/README.md` |
| `requirements.txt` | the Python environment in which the shipped outputs were produced |
| `LICENSE`, `LICENSES.md` | the licences |

Run the scripts from this folder, not from `scripts/`. No script writes outside this folder.
