#!/bin/bash
# Restart of n11v2.sh after a duplicate start (n11.sh survived TaskStop): N11a, N12
# diagnostics (5 min), then N11b and N11c. Export of 6c6f60dd, serial.
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
echo "n11v3 pid $$" >> $SPU/mx/night_status.txt
cd $SPU/mx/run
echo "N11 start $(date)" >> $SPU/mx/night_status.txt
# a. pure 2-D control plant: H-infinity Q and P, H2 o and oco, light
timeout 1800 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); [gi,h2] = mx_c1_exact(); c = {'C1_heat2d_scalar_w',struct(),'hinf','Q',gi; 'C1_heat2d_scalar_w',struct(),'hinf','P',gi; 'C1_heat2d_scalar_w',struct(),'h2','o',h2; 'C1_heat2d_scalar_w',struct(),'h2','oco',h2}; R = mx_exec2(c,'light'); save('$SP/N11a.mat','R','gi','h2');" > $SPU/mx/N11a.log 2>&1
echo "N11a done $(date)" >> $SPU/mx/night_status.txt
# N12. subsystems of S17 and S18, light
timeout 1800 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); ts = 1.034726; r = pi^2; c = {'D1_transport_1d',struct(),'stab','P',true; 'D2_s17_x2_alone',struct('r2',0),'stab','P',true; 'D2_s17_x2_alone',struct('r2',0.9*r),'stab','P',true; 'D3_s17_decoupled',struct('r2',0.9*r),'stab','P',true; 'D4_s17_heat_source',struct('r2',0.9*r),'stab','P',true; 'D4_s17_heat_source',struct('r2',1.1*r),'stab','P',false; 'D5_dde_scalar',struct('tau',0.1),'stab','P',true; 'D5_dde_scalar',struct('tau',0.5),'stab','P',true; 'D5_dde_scalar',struct('tau',0.9*ts),'stab','P',true; 'D5_dde_scalar',struct('tau',1.1*ts),'stab','P',false; 'S18_lib43_interior_delay',struct('tau',0.5,'sgn',-1,'a0',0,'a1',0.5),'stab','P',true}; R = mx_exec3(c,'light'); save('$SP/N12.mat','R');" > $SPU/mx/N12.log 2>&1
echo "N12 done $(date)" >> $SPU/mx/night_status.txt
# b. heavy tier on the other edge cascades at 0.9 of the margin
timeout 5400 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'S17_transport_edge_diffusion',struct('r2',0.9*pi^2),'stab','P',true; 'S8_two_edges',struct('r3',0.9*pi^2),'stab','P',true}; R = mx_exec2(c,'heavy'); save('$SP/N11b.mat','R');" > $SPU/mx/N11b.log 2>&1
echo "N11b done $(date)" >> $SPU/mx/night_status.txt
# c. heavy tier on the H2 coercive form, S24 (certified +12.3% at light)
timeout 3600 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'S24_hinf_ode_2d_edge',struct('a',1),'h2','oco',0.3377221541}; R = mx_exec2(c,'heavy'); save('$SP/N11c.mat','R');" > $SPU/mx/N11c.log 2>&1
echo "N11c done $(date)" >> $SPU/mx/night_status.txt
