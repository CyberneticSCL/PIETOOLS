#!/bin/bash
# cmp_ident.sh BEFORE.csv AFTER.csv
# Row-by-row identity of two harness results files on every column except
# 'seconds' (column 8) and the timings of the note (MOSEK cpu=, bisection step (Ns)), keyed by library,id,settings; rows of BEFORE absent
# from AFTER are counted only when AFTER covers that key set (M2L: a subset).
awk -F, '
NR==FNR { if (FNR>1) { k=$1","$2","$3; $8=""; t=$0; gsub(/cpu=[0-9.]+s/,"cpu=",t); gsub(/:[CFU]\([0-9]+s\)/,":()",t); a[k]=t } ; next }
FNR>1 { k=$1","$2","$3; $8=""; t=$0; gsub(/cpu=[0-9.]+s/,"cpu=",t); gsub(/:[CFU]\([0-9]+s\)/,":()",t); $0=t; n++
        if (!(k in a)) { print "NEW  " k " | " $0; nn++; next }
        if (a[k]!=$0) { d++; print "DIFF " k; print "   before: " a[k]; print "   after:  " $0 } else s++ }
END { printf "rows %d  identical %d  different %d  new %d\n", n, s, d, nn }' OFS=, "$1" "$2"
