#!/bin/bash
# P3S_FINISH.SH - 11/28/2026, Phase 3 of the _sop rename: after the batches of the tree after the
# rename (out_after). The reference is the run of 11/27 on abfb0e83 (scratchpad p3/out_new), whose
# toolbox code equals 91ea628e, the parent of the rename (only documents and records changed in
# between); two control rows on the export of 91ea628e (out_ctrl) check that the environment did
# not move. (1) SDP digests per row: reference vs out_after, reference vs out_ctrl; (2) DEMO
# workspaces; (3) F1 and F2 through p0_compare: pass 7 of C:/wt1127/new copied into the export
# after the rename as pass 8. Each MATLAB process stays far below 20 min.
S=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3s
SW='C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3s'
REF=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3/out_new
REFW='C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3/out_new'
ML="/c/Program Files/MATLAB/R2025b/bin/matlab.exe"
B=testing_codes/private/baselines/2026-11-02_opvar_phase0
echo "$(date '+%F %T') finish: start" >> "$S/drive.log"
perl "$S/p3_cmp.pl" "$REF" "$S/out_after" > "$S/cmp_sdp_all.txt" 2>&1
perl "$S/p3_cmp.pl" "$REF" "$S/out_ctrl" > "$S/cmp_sdp_ctrl.txt" 2>&1
"$ML" -batch "addpath('$SW'); p3_cmpws('$REFW/ws','$SW/out_after/ws','$SW/cmp_ws.txt')" > "$S/cmp_ws.console.txt" 2>&1
for F in F1 F2; do
  rm -rf "/c/wt1128sop/after/PIETOOLS/$B/$F/pass8"
  cp -r "/c/wt1127/new/PIETOOLS/$B/$F/pass7" "/c/wt1128sop/after/PIETOOLS/$B/$F/pass8"
done
"$ML" -sd "C:\wt1128sop\after\PIETOOLS\testing_codes\private\baselines\2026-11-02_opvar_phase0\private" -batch "setenv('P0_ROOT','C:\wt1128sop\after\PIETOOLS'); p0_compare(7,8)" > "$S/cmp_f12.console.txt" 2>&1
cp "/c/wt1128sop/after/PIETOOLS/$B/compare_pass7_pass8.txt" "$S/cmp_f12.txt" 2>/dev/null
echo "$(date '+%F %T') finish: all done" >> "$S/drive.log"
