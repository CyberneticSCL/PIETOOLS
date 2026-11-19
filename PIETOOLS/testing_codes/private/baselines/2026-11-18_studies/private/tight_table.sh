#!/bin/bash
# tight_table.sh  prints the markdown table of the B1 tightness bisections:
# one line per system and candidate T0-T3, from results_T<k>_<system>.csv
# of this folder: reach (largest certified parameter value), reach / margin,
# decision variables of the last certified program, solves, seconds, status.
# MMP, 11/18/2026.
D="$(cd "$(dirname "$0")" && pwd)"
declare -A MARGIN=( [kiss]=19.7392088 [para]=4.9348022 [s18]=1.034726 [fo2s]=2.881401 [tl]=13.686938 [age]=2.302685 [s17]=9.8696044 )
declare -A CASE=( [kiss]=34 [para]=37 [s18]=43 [fo2s]=57 [tl]=64 [age]=66 [s17]=72 )
LABEL=( "T0 light" "T1 dD 0" "T2 dD 0, Q 0 deg-raised" "T3 dD 1, Q 0 deg-raised" )
echo "| system (case) | candidate | status | reach | reach / margin | variables | solves | seconds |"
echo "|---|---|---|---|---|---|---|---|"
for s in kiss para age s18 s17 fo2s tl; do
  for k in 0 1 2 3; do
    f="$D/results_T${k}_${s}.csv"
    [ -f "$f" ] || continue
    awk -F, -v s="$s" -v c="${CASE[$s]}" -v m="${MARGIN[$s]}" -v lab="${LABEL[$k]}" 'NR>1 {
      r = ($5 ~ /^-?[0-9.eE+-]+$/) ? sprintf("%.6g", $5) : $5;
      q = ($5 ~ /^-?[0-9.eE+-]+$/) ? sprintf("%.4f", $5/m) : "";
      printf "| %s (%s) | %s | %s | %s | %s | %s | %s | %s |\n", s, c, lab, $4, r, q, $9, $7, $8 }' "$f"
  done
done
