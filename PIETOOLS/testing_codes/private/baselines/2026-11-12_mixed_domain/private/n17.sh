#!/bin/bash
# Rerun of N15a (killed by a MATLAB access violation after its first case),
# one MATLAB process per case: the systems certified at 0.9 of the margin at
# 0.98 (stable) and 1.02 (unstable), light, export of 6c6f60dd.
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
echo "N17 start pid $$ $(date)" >> $SPU/mx/night_status.txt
cd $SPU/mx/run
K="k0 = fzero(@(k) k*tan(k)-1,[0.1 1.5]); r4 = pi^2+k0^2; a26 = pi^2/4-0.5;"
CASES=(
"{'S6_indomain_1d_2d',struct('g',0.98*pi^2),'stab','P',true}"
"{'S7_ode_2d_indomain',struct('a',1,'b',1,'g',0.98),'stab','P',true}"
"{'S9_ode_2d_flux',struct('a',0.98,'c',1),'stab','P',true}"
"{'S26_lib42_ode_heat_delay',struct('a',0.98*a26),'stab','P',true}"
"{'S4_wentzell',struct('c',1,'r',1.02*r4),'stab','P',false}"
"{'S6_indomain_1d_2d',struct('g',1.02*pi^2),'stab','P',false}"
"{'S7_ode_2d_indomain',struct('a',1,'b',1,'g',1.02),'stab','P',false}"
"{'S9_ode_2d_flux',struct('a',1.02,'c',1),'stab','P',false}"
"{'S26_lib42_ode_heat_delay',struct('a',1.02*a26),'stab','P',false}"
"{'S1_ode_heat1d_bc',struct('a',0.02),'stab','P',true}"
"{'S2_ode_heat2d_edge',struct('a',0.02),'stab','P',true}"
"{'S5_chain_ode_1d_2d',struct('a',0.02),'stab','P',true}"
"{'S1_ode_heat1d_bc',struct('a',-0.02),'stab','P',false}"
"{'S2_ode_heat2d_edge',struct('a',-0.02),'stab','P',false}"
"{'S5_chain_ode_1d_2d',struct('a',-0.02),'stab','P',false}"
)
: > $SPU/mx/N17.log
i=0
for c in "${CASES[@]}"; do
    i=$((i+1))
    echo "=== N17 case $i: $c" >> $SPU/mx/N17.log
    timeout 900 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); $K R = mx_exec3($c,'light'); save('$SP/N17_$i.mat','R');" >> $SPU/mx/N17.log 2>&1
    echo "=== N17 case $i exit $?" >> $SPU/mx/N17.log
done
echo "N17 done $(date)" >> $SPU/mx/night_status.txt
