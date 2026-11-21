#!/bin/bash
# review.sh GROUP_TSV - coordinator check of one group: for each condense file, the TAG of its
# note, cpass.pl check, codeident against HEAD, line-ending state, whitespace-only changes; for
# asnew files, no dated note left but the Initial coding line; review/deferred files untouched.
T=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/cpass1121
RC=cleaning/private/release_comment_pass
bad=0
while IFS=$'\t' read -r f cls rest; do
  case "$cls" in
  condense)
    tag=$(tr -d '\r' < "$f" | sed -n 's/^% \(.*\): Changes since PIETOOLS 2025: .*/\1/p' | head -1)
    if [ -z "$tag" ]; then echo "NOTE MISSING  $f"; bad=1; continue; fi
    r=$(perl $T/cpass.pl check "$f" "$tag" | tr '\n' ' ')
    case "$r" in PASS*) ;; *) echo "$r"; bad=1;; esac ;;
  asnew)
    n=$(tr -d '\r' < "$f" | grep -c -E '%.*\b[A-Z]{2,3}( *\(CC\))?( *, *[A-Z]{2,3})*,? *[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}|\(was\)')
    i=$(tr -d '\r' < "$f" | grep -c -i 'initial coding')
    [ "$n" -gt "$i" ] && { echo "ASNEW: $n dated lines (Initial coding $i)  $f"; bad=1; } ;;
  review|deferred)
    git diff --quiet -- "$f" || { echo "TOUCHED ($cls)  $f"; bad=1; }
    continue ;;
  esac
  c=$(bash $RC/codeident.sh "$f" | head -1); case "$c" in OK*) ;; *) echo "CODE $c"; bad=1;; esac
  e=$(git ls-files --eol -- "$f" | awk '{print $2}'); [ "$e" = "w/crlf" ] || echo "EOL $e  $f"
  a=$(git diff --numstat -- "$f"); b=$(git diff -w --numstat -- "$f"); [ "$a" = "$b" ] || { echo "WHITESPACE  $f"; bad=1; }
done < "$1"
[ $bad -eq 0 ] && echo "review $(basename $1): all PASS"
