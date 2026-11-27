#!/bin/bash
# REN_CHECK.SH BASE NEW PATHS DATE OUT - checks 3a-3c of Phase 1 (PLAN.md), 11/26/2026, written
# apart from ren_apply.pl: the map is a sed script generated from rename_table.tsv.
#   BASE, NEW: the PIETOOLS folder before and after; PATHS: "old<TAB>new" moves (ren_apply.pl
#   paths); DATE: mm/dd/yyyy of the D1 lines; OUT: folder for the expected files and reports.
#   0. every file outside the moves is byte-identical, and each moved file exists;
#   3b. code identity: codeident(NEW file) = codeident(sed map of the BASE file);
#   3c. comment identity (ren_cmp.pl): NEW = sed map of BASE line by line, except the D1 lines,
#       the records (header-log and (was) lines, kept as in BASE), blank runs, and an article
#       that ends a line before a mapped cop... token;
#   3a. no token naming a container outside the records and the retired names (ren_cmp.pl).
B=$1; N=$2; PATHS=$3; DATE=$4; OUT=$5
HERE=$(dirname "$0"); TAB=$HERE/../rename_table.tsv
CI=$HERE/../../release_comment_pass/codeident.awk
mkdir -p "$OUT/expect"; rm -f "$OUT/pairs.tsv"
[ -f "$OUT/changed.txt" ] || { echo "OUT/changed.txt (ren_apply.pl report) missing"; exit 2; }
# the sed map: the article first, then each name of the table as a whole word, in lower,
# upper and capitalized form; (retire) rows stay
awk -F'\t' 'NR>1 && $6 !~ /retire/ { o=tolower($1); n=tolower($6);
    printf "s/\\b%s\\b/%s/g\n", o, n; printf "s/\\b%s\\b/%s/g\n", toupper(o), toupper(n);
    printf "s/\\b%s\\b/%s/g\n", toupper(substr(o,1,1)) substr(o,2), toupper(substr(n,1,1)) substr(n,2);
    if (o ~ /^cop/ && n ~ /^op/) art = art (art ? "\\|" : "") o "\\|" toupper(o) "\\|" toupper(substr(o,1,1)) substr(o,2) }
  END { printf "s/\\b\\([Aa]\\)\\([[:space:]]\\+[\x27\"`(]\\?\\)\\(%s\\)\\b/\\1n\\2\\3/g\n", art }' "$TAB" > "$OUT/map_body.sed"
{ tail -1 "$OUT/map_body.sed"; head -n -1 "$OUT/map_body.sed"; } > "$OUT/map.sed"
echo "sed map: $(wc -l < "$OUT/map.sed") rules"
# 0. files outside the moves
declare -A NEWOF; while IFS=$'\t' read -r o n; do NEWOF[$o]=$n; done < "$PATHS"
nb=0; nm=0
while IFS= read -r f; do
  n=${NEWOF[$f]:-$f}
  if [ ! -f "$N/$n" ]; then echo "MISSING $n"; nm=$((nm+1)); continue; fi
  if [ "$n" == "$f" ] && ! grep -qxF "$f" "$OUT/changed.txt"; then cmp -s "$B/$f" "$N/$n" || { echo "CHANGED-UNLISTED $f"; nb=$((nb+1)); }; fi
done < <(cd "$B" && find . -type f | sed 's#^\./##' | grep -vE '^(cleaning|orphaned_file|SOSTOOLS400)/')
echo "0. unlisted files changed: $nb, missing: $nm"
# 3b. code identity on the changed files
ok=0; bad=0
while IFS= read -r n; do
  o=$n; for k in "${!NEWOF[@]}"; do [ "${NEWOF[$k]}" == "$n" ] && o=$k && break; done
  mkdir -p "$(dirname "$OUT/expect/$n")"
  sed -f "$OUT/map.sed" "$B/$o" > "$OUT/expect/$n"
  if [ "$(awk -f "$CI" "$OUT/expect/$n")" == "$(awk -f "$CI" "$N/$n")" ]; then ok=$((ok+1)); else bad=$((bad+1)); echo "3b DIFF $n"; fi
  printf '%s\t%s\n' "$o" "$n" >> "$OUT/pairs.tsv"
done < "$OUT/changed.txt"
echo "3b. code identity: $ok OK, $bad DIFF"
perl "$HERE/ren_cmp.pl" "$B" "$N" "$OUT" "$DATE"
