#!/bin/bash
# After N14: (a) tightness and soundness near each margin, light: the
# systems certified at 0.9 run at 0.98 (stable) and 1.02 (unstable);
# (b) heavy, the tier that certifies S3, at 1.02 and 1.1 r* (unstable).
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
echo "n15 pid $$ waiting $(date)" >> $SPU/mx/night_status.txt
until grep -q "N14 done" $SPU/mx/night_status.txt 2>/dev/null; do sleep 30; done
cd $SPU/mx/run
timeout 2400 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); k0 = fzero(@(k) k*tan(k)-1,[0.1 1.5]); r4 = pi^2+k0^2; a26 = pi^2/4-0.5; c = {}; for f=[0.98 1.02], st = f<1; c(end+1,:) = {'S4_wentzell',struct('c',1,'r',f*r4),'stab','P',st}; c(end+1,:) = {'S6_indomain_1d_2d',struct('g',f*pi^2),'stab','P',st}; c(end+1,:) = {'S7_ode_2d_indomain',struct('a',1,'b',1,'g',f),'stab','P',st}; c(end+1,:) = {'S9_ode_2d_flux',struct('a',f,'c',1),'stab','P',st}; c(end+1,:) = {'S26_lib42_ode_heat_delay',struct('a',f*a26),'stab','P',st}; end; for a=[0.02 -0.02], c(end+1,:) = {'S1_ode_heat1d_bc',struct('a',a),'stab','P',a>0}; c(end+1,:) = {'S2_ode_heat2d_edge',struct('a',a),'stab','P',a>0}; c(end+1,:) = {'S5_chain_ode_1d_2d',struct('a',a),'stab','P',a>0}; end; R = mx_exec3(c,'light'); save('$SP/N15a.mat','R');" > $SPU/mx/N15a.log 2>&1
echo "N15a done $(date)" >> $SPU/mx/night_status.txt
timeout 2400 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'S3_edge1d_heat2d',struct('r',1.02*pi^2),'stab','P',false; 'S3_edge1d_heat2d',struct('r',1.1*pi^2),'stab','P',false}; R = mx_exec3(c,'heavy'); save('$SP/N15b.mat','R');" > $SPU/mx/N15b.log 2>&1
echo "N15b done $(date)" >> $SPU/mx/night_status.txt
