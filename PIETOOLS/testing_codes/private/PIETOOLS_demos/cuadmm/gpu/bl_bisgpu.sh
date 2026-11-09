#!/bin/bash
# bl_bisgpu.sh -- cost of ONE cuADMM bisection step.
#
# The step is not "solve to tolerance": measured on h2c_rd1, neither the feasible
# nor the infeasible pinned program converges, because cuADMM's stopping rule
# needs the duality GAP and the gap does not close on these. What separates them
# is the primal residual PLATEAU -- the infeasible program sticks at
# pinf ~ 1.5e-02 from about iteration 500 and never improves, while the feasible
# one keeps descending (9.7e-3 -> 2.1e-3 -> 5.1e-4 -> 1.9e-5). So a bisection
# step is a FIXED ITERATION BUDGET followed by a threshold on pinf.
#
# 5000 iterations separates the two by 30x on that case. Reported here: pinf at
# 1000/2000/5000 so the budget can be chosen from data rather than asserted, and
# s/iteration so the cost is comparable across machines and load.
set -u
EXE=${CUADMM_EXE:-/home/mpeet/solvers/cuADMM/build/cuadmm_exe}   # override with CUADMM_EXE
B="$1"; OUT="$2"; BUDGET=${3:-5000}
printf 'lab\tbudget\tsolve_s\ts_per_it\tpinf_1000\tpinf_2000\tpinf_final\tdinf_final\tgap_final\tdobj_final\n' > "$OUT"
for d in "$B"/*/; do
  lab=$(basename "$d"); [ -f "$d/blk.txt" ] || continue
  log="$d/step_${BUDGET}.log"
  LD_LIBRARY_PATH=/usr/lib/wsl/lib "$EXE" "$d" 1e-12 "$BUDGET" 15 0 > "$log" 2>&1
  sv=$(grep -oP '(?<= solve )[0-9.eE+-]+' "$log" | tail -1)
  read p1 p2 pf df gf dof <<<$(awk '/^ *[0-9]+ \|/{gsub(/\|/," ");
      if($1==1000)p1=$2; if($1==2000)p2=$2; pf=$2; df=$3; gf=$6; dof=$5}
      END{print (p1==""?"-":p1), (p2==""?"-":p2), pf, df, gf, dof}' "$log")
  # "$log" was missing here: awk read empty stdin and every residual column came
  # out blank with no error. results/2026-09-24/bl_bisgpu.tsv shows it; the values
  # the report cites were read by hand from results/2026-09-24/logs/bisect/.
  spi=$(awk -v s="${sv:-0}" -v b="$BUDGET" 'BEGIN{if(b>0)printf "%.5f", s/b; else print "-"}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lab" "$BUDGET" "${sv:-}" "$spi" "$p1" "$p2" "$pf" "$df" "$gf" "$dof" >> "$OUT"
  echo "BG $lab budget=$BUDGET solve=${sv:-NA}s s/it=$spi pinf=$pf gap=$gf"
done
echo "BGDONE"
