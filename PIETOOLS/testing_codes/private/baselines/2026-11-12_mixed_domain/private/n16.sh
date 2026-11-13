#!/bin/bash
# After N15: S8 at heavy, 0.9 r3*, only if N11b (cap 90 min, spent mostly on
# S17 at dw 2) left no S8 result. Export of 6c6f60dd.
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
echo "n16 pid $$ waiting $(date)" >> $SPU/mx/night_status.txt
until grep -q "N15b done" $SPU/mx/night_status.txt 2>/dev/null; do sleep 30; done
if grep -q "^S8_two_edges" $SPU/mx/N11b.log; then
    echo "N16 skipped (S8 in N11b) $(date)" >> $SPU/mx/night_status.txt
    exit 0
fi
cd $SPU/mx/run
timeout 3600 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'S8_two_edges',struct('r3',0.9*pi^2),'stab','P',true}; R = mx_exec3(c,'heavy'); save('$SP/N16.mat','R');" > $SPU/mx/N16.log 2>&1
echo "N16 done $(date)" >> $SPU/mx/night_status.txt
