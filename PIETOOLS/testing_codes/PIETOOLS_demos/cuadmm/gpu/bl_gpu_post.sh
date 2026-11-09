#!/bin/bash
# bl_gpu_post.sh -- last-iteration residuals for EVERY cuADMM run, converged or
# capped. The summary block only exists on a converged run, so a capped case
# would otherwise contribute nothing but "capped" -- yet the capped runs carry
# the most informative number in this whole arm: whether the run stalled on
# FEASIBILITY or on the duality GAP.
#
# cuADMM stops on max(pinf,dinf,relgap) < tol. Measured on h2oco_rd1 (m=56):
# pinf 1.36e-05 and dinf 1.42e-05 are both already inside 1e-4 while relgap sits
# at 4.56e-04 and falls by ~1e-8 per 100 iterations. So "did not converge" means
# "did not close the gap", not "did not reach feasibility" -- a distinction that
# decides whether bisection (a feasibility form) rescues the objective classes.
#
# Trajectory row: iter | pinf dinf | pobj dobj gap | time | sigma
DMP="$1"; OUT="$2"
printf 'id\ttol\titers\tlast_pinf\tlast_dinf\tlast_gap\tlast_pobj\tlast_dobj\tlast_t\tlimiter\n' > "$OUT"
for d in "$DMP"/*/; do
  id=$(basename "$d")
  for tol in 1e-4 1e-6; do
    log="$d/cu_${tol}.log"; [ -f "$log" ] || continue
    awk -v id="$id" -v tol="$tol" '
      /^ *[0-9]+ \|/ { gsub(/\|/," "); it=$1; pi=$2; di=$3; po=$4; dobj=$5; gp=$6; tt=$7 }
      END{
        if (it=="") { exit }
        lim = "gap"
        if (pi+0 >= gp+0 && pi+0 >= di+0) lim="pinf"
        else if (di+0 >= gp+0) lim="dinf"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", id,tol,it,pi,di,gp,po,dobj,tt,lim
      }' "$log" >> "$OUT"
  done
done
echo "POSTDONE $(wc -l < "$OUT") rows"
