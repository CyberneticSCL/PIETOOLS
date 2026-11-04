# lb_augment.awk, 10/28/2026: comparison.csv (lb_compare.awk) plus, per
# path, the construction/solve split and equality count (lb_split.awk) and
# a second verdict that certifies by MOSEK's solution status.
# Usage: awk -F, -f lb_augment.awk split_stock.csv split_sop.csv \
#            rows_stock.csv rows_sop.csv comparison.csv > comparison_full.csv
# Added columns, per path (stock_, sop_): sdps, solve_seconds,
# construct_seconds, equalities (as in lb_split.awk), mosek (solution
# status of the row's last SDP); then verdict_mosek: the verdict of
# lb_compare.awk with "certified" read as "the last SDP's solution status
# is OPTIMAL", for rows with one executive call on both paths; bisection
# rows keep the feasratio verdict (each step was classified by it).
# Summary lines prefixed '#' on stderr.
function cert(s) { return s=="OPTIMAL" }
FILENAME ~ /split_stock/ { if (FNR>1) { k=$1"|"$2"|"$3; a1[k]=$7","$9","$10","$11 } ; next }
FILENAME ~ /split_sop/   { if (FNR>1) { k=$1"|"$2"|"$3; b1[k]=$7","$9","$10","$11 } ; next }
FILENAME ~ /rows_stock/  { i=index($1,":"); k=substr($1,1,i-1)"|"substr($1,i+1)"|"$2; a2[k]=$7; next }
FILENAME ~ /rows_sop/    { i=index($1,":"); k=substr($1,1,i-1)"|"substr($1,i+1)"|"$2; b2[k]=$7; next }
FNR==1 { print $0",stock_sdps,stock_solve_seconds,stock_construct_seconds,stock_equalities,sop_sdps,sop_solve_seconds,sop_construct_seconds,sop_equalities,stock_mosek,sop_mosek,verdict_mosek"; next }
{
  k=$1"|"$2"|"$3; v=$16
  sa=(k in a1 ? a1[k] : "0,NaN,NaN,NaN"); sb=(k in b1 ? b1[k] : "0,NaN,NaN,NaN")
  ma=(k in a2 ? a2[k] : "-"); mb=(k in b2 ? b2[k] : "-")
  if (v ~ /error|skipped/) w=v
  else if ($5=="yes") w=v
  else {
    A=cert(ma); B=cert(mb)
    if (A && !B) w="stock only"; else if (!A && B) w="sop only"; else if (!A && !B) w="neither"
    else {
      x=$13+0; y=$9+0; rel=(x-y)/((y<0?-y:y)+1e-12)
      if ($4 ~ /^stability/) w="same"; else if (rel<1e-6 && rel>-1e-6) w="same"
      else w=(x<y ? "sop tighter" : "sop looser")
    }
  }
  print $0","sa","sb","ma","mb","w
  n1[$3" "v]++; n2[$3" "w]++
}
END {
  for (q in n1) print "# feasratio rule  " q ": " n1[q] > "/dev/stderr"
  for (q in n2) print "# MOSEK rule      " q ": " n2[q] > "/dev/stderr"
}
