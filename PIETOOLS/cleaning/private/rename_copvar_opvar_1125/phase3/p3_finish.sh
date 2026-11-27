#!/bin/bash
# p3_finish.sh - 11/27/2026: after p3_drive.sh. (1) SDP digests per row, old vs new; (2) DEMO
# workspaces; (3) F1 and F2 through p0_compare: the old pass 7 copied into the new export as pass 8;
# (4) batch X (3-D heat LPI) on both trees, then its digests; (5) class-operation timings,
# alternating the trees twice. Each MATLAB process stays far below 20 min.
S=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3
SW='C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3'
ML="/c/Program Files/MATLAB/R2025b/bin/matlab.exe"
B=testing_codes/private/baselines/2026-11-02_opvar_phase0
echo "$(date '+%F %T') finish: start" >> "$S/drive.log"
perl "$S/p3_cmp.pl" "$S/out_old" "$S/out_new" > "$S/cmp_sdp.txt" 2>&1
"$ML" -batch "addpath('$SW'); p3_cmpws('$SW/out_old/ws','$SW/out_new/ws','$SW/cmp_ws.txt')" > "$S/cmp_ws.console.txt" 2>&1
for F in F1 F2; do
  rm -rf "/c/wt1127/new/PIETOOLS/$B/$F/pass8"
  cp -r "/c/wt1127/old/PIETOOLS/$B/$F/pass7" "/c/wt1127/new/PIETOOLS/$B/$F/pass8"
done
"$ML" -sd "C:\\wt1127\\new\\PIETOOLS\\testing_codes\\private\\baselines\\2026-11-02_opvar_phase0\\private" -batch "setenv('P0_ROOT','C:\\wt1127\\new\\PIETOOLS'); p0_compare(7,8)" > "$S/cmp_f12.console.txt" 2>&1
cp "/c/wt1127/new/PIETOOLS/$B/compare_pass7_pass8.txt" "$S/cmp_f12.txt" 2>/dev/null
echo "$(date '+%F %T') finish: comparisons done" >> "$S/drive.log"
bash "$S/p3_drive.sh" "old new" "X"
perl "$S/p3_cmp.pl" "$S/out_old" "$S/out_new" > "$S/cmp_sdp_all.txt" 2>&1
for t in old new old new; do
  "$ML" -batch "addpath('$SW'); t3_classops('C:/wt1127/$t')" 2>&1 | grep -E "^tree|median|decision" >> "$S/classops.txt"
done
echo "$(date '+%F %T') finish: all done" >> "$S/drive.log"
