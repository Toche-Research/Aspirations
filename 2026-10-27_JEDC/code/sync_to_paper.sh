#!/bin/sh
# Copy the generated tables and figures from this repository into the paper's own directories,
# which are not published. Run from code/ after regenerating:
#
#     sh sync_to_paper.sh
#
# The scripts write into out/tables/ and out/figures/, so that everything a replicator produces
# stays inside the repository. This step is only for building the PDF of the paper.
#
# A group that has nothing to copy reports zero rather than aborting: not every paper in this
# project generates both tables and figures, and a cleared out/ is a legitimate state.
set -eu
cd "$(dirname "$0")"

TABDEST=../tex/tables
FIGDEST=../tex/figures
MACDEST=../tex/macros
TIKZDEST=../tex/tikz
mkdir -p "$TABDEST" "$FIGDEST" "$MACDEST" "$TIKZDEST"

copy_group() {          # $1 source directory, $2 extension, $3 destination, $4 what they are
    n=0
    for f in "$1"/*"$2"; do
        [ -e "$f" ] || break
        cp "$f" "$3/"
        n=$((n + 1))
    done
    echo "copied $n $4 to $3/"
}

copy_group out/tables .tex "$TABDEST" "table bodies"
copy_group out/figures .pdf "$FIGDEST" "figures"
copy_group out/macros .tex "$MACDEST" "macro files"
copy_group out/tikz .tex "$TIKZDEST" "TikZ fragments"
