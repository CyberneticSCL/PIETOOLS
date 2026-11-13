#!/bin/bash
# After N13c: the default H2 form 'c' (and 'cco') on the control C1, which
# has no K, so 'c' is posed as itself. Light, export of 6c6f60dd.
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
echo "n14 pid $$ waiting $(date)" >> $SPU/mx/night_status.txt
until grep -q "N13c done" $SPU/mx/night_status.txt 2>/dev/null; do sleep 30; done
cd $SPU/mx/run
timeout 1800 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); [gi,h2] = mx_c1_exact(); c = {'C1_heat2d_scalar_w',struct(),'h2','c',h2; 'C1_heat2d_scalar_w',struct(),'h2','cco',h2}; R = mx_exec2(c,'light'); save('$SP/N14.mat','R');" > $SPU/mx/N14.log 2>&1
echo "N14 done $(date)" >> $SPU/mx/night_status.txt
