#!/bin/bash
# Phase 0 driver. Usage: bash run_phase0.sh now|after0430|rest|c81now|c81after0430
# Tree eb5d0600 (export p0_eb5d0600), passes 1 and 2:
#   now        F1 and F2, pass 1 (estimated 10-14 min; inside the 20-minute rule)
#   after0430  waits until 04:30, then F3 pass 1 (light, heavy), F4 (DEMO1-11),
#              F1 and F2 pass 2, F3 pass 2 (light), and the pass comparison
#              (estimated 75 min, worst case set by the caps about 4.4 h)
#   rest       11/02 06:30: light rows after PDE 50 in both passes, PDE 50 skipped
#              (2,421 s per solve), then the comparison (estimated 9-10 min)
# Tree c81c3d04 (export p0_c81c3d04), passes 3 and 4:
#   c81now        F1 and F2, pass 3, compared with pass 1 (estimated 8-12 min;
#                 the caps total 20 min)
#   c81after0430  waits until 04:30, then F3 light passes 3 and 4, F4 pass 3, and
#                 the comparisons 1-3 and 3-4 (estimated 38 min: 16 min per F3
#                 pass, 5 min F4; worst case set by the caps about 3.3 h)
# Each step is one MATLAB process with a time cap. Create the file STOP in this
# folder to skip the steps that have not started. Log: runs_phase0.log.
# MMP, 11/05/2026: DEMO11 is DEMO11_2DControl, which replaced
#                  DEMO11_CdopvarControl in 36f50e26 (same LPI; B2 -> Bu and
#                  the variables from PIE.space). The F4 comparison pairs
#                  passes by file name, so passes 1-4 have no partner for it.
P0W="C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/baselines/2026-11-02_opvar_phase0"
LOG="$P0W/runs_phase0.log"
NOWAIT="$2"
mkdir -p "$P0W/logs"
# DEMO11 renamed; the third line of DEMOS was                               # MMP, 11/05/2026
#      DEMO7_observer_based_control DEMO8_H2norm DEMO9_L2_gain_analysis DEMO10_PBC DEMO11_CdopvarControl" # MMP, 11/05/2026 (was)
DEMOS="DEMO1_stability_and_simulation DEMO2_volterra_operator_norm DEMO3_poincare_inequality
       DEMO4_stability_parameter_bisection DEMO5_Hinf_optimal_estimator DEMO6_Hinf_optimal_control
       DEMO7_observer_based_control DEMO8_H2norm DEMO9_L2_gain_analysis DEMO10_PBC DEMO11_2DControl"
step() {   # name, cap in seconds, MATLAB command
  if [ -f "$P0W/STOP" ]; then echo "$(date '+%F %T') STOP present: $1 not run" >> "$LOG"; return; fi
  echo "$(date '+%F %T') start $1 (cap $2 s, P0_ROOT $P0_ROOT): $3" >> "$LOG"
  timeout "$2" matlab -sd "$P0W/private" -batch "$3" > "$P0W/logs/$1.txt" 2>&1
  rc=$?
  echo "$(date '+%F %T') end $1 exit $rc" >> "$LOG"
}
wait0430() {   # a second argument 'nowait' starts at once (11/04 02:05, at the maintainer's request)
  if [ "$NOWAIT" = "nowait" ]; then echo "$(date '+%F %T') starting without waiting for 04:30" >> "$LOG"; return; fi
  now=$(date +%s); target=$(date -d "04:30" +%s)
  if [ "$now" -lt "$target" ]; then
    echo "$(date '+%F %T') waiting until 04:30 ($((target-now)) s)" >> "$LOG"
    sleep $((target-now))
  fi
}
case "$1" in
  now|after0430|rest) export P0_ROOT='C:\Users\mpeet\wt\p0_eb5d0600' ;;
  c81now|c81after0430) export P0_ROOT='C:\Users\mpeet\wt\p0_c81c3d04' ;;
esac
case "$1" in
  now)
    step F1_pass1 1200 "p0_f1(1)"
    step F2_pass1 1200 "p0_f2(1)"
    ;;
  after0430)
    wait0430
    step F3_pass1 4500 "p0_f3(1,{'light','heavy'})"
    for d in $DEMOS; do
      step "F4_pass1_$d" 600 "p0_f4one('$d',1)"
    done
    step F1_pass2 1200 "p0_f1(2)"
    step F2_pass2 1200 "p0_f2(2)"
    step F3_pass2 1800 "p0_f3(2,{'light'})"
    step compare 600 "p0_compare(1,2)"
    echo "$(date '+%F %T') after0430 finished" >> "$LOG"
    ;;
  rest)
    step F3_light_rest 1200 "p0_f3(1,{'light'}); p0_f3(2,{'light'}); p0_compare(1,2)"
    ;;
  c81now)
    step F1_pass3 690 "p0_f1(3)"
    step F2_pass3 390 "p0_f2(3)"
    step compare13_F12 120 "p0_compare(1,3)"
    echo "$(date '+%F %T') c81now finished" >> "$LOG"
    ;;
  c81after0430)
    wait0430
    step F3_pass3 1800 "p0_f3(3,{'light'})"
    step F3_pass4 1800 "p0_f3(4,{'light'})"
    for d in $DEMOS; do
      step "F4_pass3_$d" 600 "p0_f4one('$d',3)"
    done
    step compare13 600 "p0_compare(1,3)"
    step compare34 600 "p0_compare(3,4)"
    echo "$(date '+%F %T') c81after0430 finished" >> "$LOG"
    ;;
  *) echo "usage: bash run_phase0.sh now|after0430|rest|c81now|c81after0430"; exit 1 ;;
esac
