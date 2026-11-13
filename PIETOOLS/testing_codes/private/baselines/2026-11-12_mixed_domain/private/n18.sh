#!/bin/bash
# N17 cases without a result: S26 (absent from MX_SYS) at 0.98 and 1.02 of
# its margin, and S7 at 1.02 (access violation in N17), through the test
# catalogue, one MATLAB process per case, light, export of 6c6f60dd.
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
echo "N18 start pid $$ $(date)" >> $SPU/mx/night_status.txt
cd $SPU/mx/run
CASES=(
"{'S26_lib42_ode_heat_delay',struct('a',0.98*(pi^2/4-0.5)),'stab','P',true}"
"{'S26_lib42_ode_heat_delay',struct('a',1.02*(pi^2/4-0.5)),'stab','P',false}"
"{'S7_ode_2d_indomain',struct('a',1,'b',1,'g',1.02),'stab','P',false}"
)
: > $SPU/mx/N18.log
i=0
for c in "${CASES[@]}"; do
    i=$((i+1))
    echo "=== N18 case $i: $c" >> $SPU/mx/N18.log
    timeout 900 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); R = mx_exec4($c,'light'); save('$SP/N18_$i.mat','R');" >> $SPU/mx/N18.log 2>&1
    echo "=== N18 case $i exit $?" >> $SPU/mx/N18.log
done
echo "N18 done $(date)" >> $SPU/mx/night_status.txt
