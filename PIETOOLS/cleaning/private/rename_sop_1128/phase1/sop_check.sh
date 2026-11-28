#!/bin/bash
# SOP_CHECK.SH BASE NEW PATHS DATE OUT - checks 0 and 3a-3c of Phase 1 (../PLAN.md), 11/28/2026,
# written apart from sop_apply.pl (adapted from ../../rename_copvar_opvar_1125/phase1/ren_check.sh):
# the map is a sed script generated from ../rename_table.tsv.
#   BASE, NEW: the PIETOOLS folder before and after; PATHS: "old<TAB>new" moves (sop_apply.pl
#   paths); DATE: mm/dd/yyyy of the D1 lines; OUT: the report folder of sop_apply.pl apply
#   (changed.txt), which receives the expected files and the reports.
#   0. every file outside the moves and outside changed.txt is byte-identical; each moved file
#      exists;
#   3b. code identity: codeident(NEW file) = codeident(sed map of the BASE file);
#   3c. comment identity (sop_cmp.pl): NEW = sed map of BASE line by line, except the D1 lines,
#       the records (header-log and (was) lines, kept as in BASE) and runs of blanks;
#   3a. no old name of the table outside the records (sop_cmp.pl).
B=$1; N=$2; PATHS=$3; DATE=$4; OUT=$5
HERE=$(dirname "$0"); TAB=$HERE/../rename_table.tsv
CI=$HERE/../../release_comment_pass/codeident.awk
mkdir -p "$OUT/expect"; rm -f "$OUT/pairs.tsv"
[ -f "$OUT/changed.txt" ] || { echo "OUT/changed.txt (sop_apply.pl report) missing"; exit 2; }
# the sed map: two text rules, then each renamed name of the table and each label of the
# executives as a whole word, in lower, upper and capitalized form
{ echo "s/PIETOOLS_\\*_sop\\b/PIETOOLS_*/g"
  echo "/^[[:space:]]*%/!s/[[:space:]]\\+'_sop'\\]/]/g"
  { awk -F'\t' 'NR>1 && $4 != $1 { print $1 "\t" $4 }' "$TAB"
    for l in stability Hinf_gain H2_norm well_posedness H2_control H2_estimator Hinf_OFcontrol \
             Hinf_control Hinf_estimator; do printf '%s_sop\t%s\n' "$l" "$l"; done
    printf '_Hinf_gain_sop\t_Hinf_gain\n_H2_norm_sop\t_H2_norm\n'
  } | awk -F'\t' '{ o=$1; n=$2;
      printf "s/\\b%s\\b/%s/g\n", o, n; printf "s/\\b%s\\b/%s/g\n", toupper(o), toupper(n);
      printf "s/\\b%s\\b/%s/g\n", toupper(substr(o,1,1)) substr(o,2), toupper(substr(n,1,1)) substr(n,2) }'
} > "$OUT/map.sed"
echo "sed map: $(wc -l < "$OUT/map.sed") rules"
# 0. files outside the moves
declare -A NEWOF; while IFS=$'\t' read -r o n; do NEWOF[$o]=$n; done < "$PATHS"
declare -A CHG; while IFS= read -r f; do CHG[$f]=1; done < "$OUT/changed.txt"
nb=0; nm=0
while IFS= read -r f; do
  n=${NEWOF[$f]:-$f}
  if [ ! -f "$N/$n" ]; then echo "MISSING $n"; nm=$((nm+1)); continue; fi
  if [ -z "${CHG[$n]}" ]; then cmp -s "$B/$f" "$N/$n" || { echo "CHANGED-UNLISTED $f"; nb=$((nb+1)); }; fi
done < <(cd "$B" && find . -type f | sed 's#^\./##' | grep -vE '^(cleaning|orphaned_file|SOSTOOLS400)/')
echo "0. unlisted files changed: $nb, missing: $nm"
# 3b. code identity on the changed files
declare -A OLDOF; for k in "${!NEWOF[@]}"; do OLDOF[${NEWOF[$k]}]=$k; done
ok=0; bad=0
while IFS= read -r n; do
  o=${OLDOF[$n]:-$n}
  mkdir -p "$(dirname "$OUT/expect/$n")"
  sed -f "$OUT/map.sed" "$B/$o" > "$OUT/expect/$n"
  if [ "$(awk -f "$CI" "$OUT/expect/$n")" == "$(awk -f "$CI" "$N/$n")" ]; then ok=$((ok+1)); else bad=$((bad+1)); echo "3b DIFF $n"; fi
  printf '%s\t%s\n' "$o" "$n" >> "$OUT/pairs.tsv"
done < "$OUT/changed.txt"
echo "3b. code identity: $ok OK, $bad DIFF"
perl "$HERE/sop_cmp.pl" "$B" "$N" "$OUT" "$DATE" "$TAB"
