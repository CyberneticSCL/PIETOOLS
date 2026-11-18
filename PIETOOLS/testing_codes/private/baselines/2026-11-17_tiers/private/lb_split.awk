# lb_split.awk, 10/28/2026: the construction/solve split and the equality
# count of every row of an lb_run results file, from the MOSEK runs in its log.
# Usage: awk -v csv=1 -f lb_mosekrows.awk LOG > rows.csv
#        awk -F, -v path=stock -f lb_split.awk rows.csv RESULTS.csv
# Output columns:
#   library,id,settings,path,status,n_solves,n_sdp,seconds,solve_seconds,
#   construct_seconds,equalities,decision_variables
# solve_seconds: MOSEK optimizer time (its log, "Optimizer - time"), summed
#   over every SDP of the row (bisection steps, slack-raise passes);
# construct_seconds: seconds - solve_seconds, i.e. building the LPI,
#   converting it to solver input, and reading the solution back;
# equalities: MOSEK "Constraints" of the row's last SDP (rows of the SeDuMi
#   equality system); NaN when the log holds no SDP for the row.
FNR==NR { k=$1"|"$2; ns[k]=$3; ts[k]=$4; mc[k]=$5; next }
FNR==1  { print "library,id,settings,path,status,n_solves,n_sdp,seconds,solve_seconds,construct_seconds,equalities,decision_variables"; next }
{
  k=$1":"$2"|"$3
  if (k in ns) { n=ns[k]; s=ts[k]; c=$8-s; m=(mc[k]=="-" ? "NaN" : mc[k]) } else { n=0; s="NaN"; c="NaN"; m="NaN" }
  printf "%s,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s,%s\n", $1, $2, $3, path, $4, $7, n, $8, s, (c=="NaN" ? c : sprintf("%.2f",c)), m, $9
}
