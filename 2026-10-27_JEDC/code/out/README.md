<!-- Describes what belongs in this directory. -->

# code/out

Everything the code produces. **This whole tree is generated** and is reproducible by rerunning the
scripts; nothing in it should be edited by hand.

| subdirectory | what it holds |
|---|---|
| `figures/` | figure PDFs, copied to `tex/figures/` |
| `tables/` | table bodies, copied to `tex/tables/` |
| `data/` | intermediate series and solution objects |

`code/sync_to_paper.sh` copies `tables/`, `figures/` and `data/` into `tex/`. It is the only thing that should
write into `tex/figures/`, `tex/tables/` or `tex/macros/`.

The outputs are committed deliberately, so the paper compiles without running anything first.
