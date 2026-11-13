#!/bin/bash
# After N11c: (a) the control C1 through the test driver (convert, ops,
# exec on C1; the whole exact layer); (b) S18 at other a0, a1 (light);
# (c) C1 at heavy. Export of 6c6f60dd, serial.
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
echo "n13 pid $$ waiting $(date)" >> $SPU/mx/night_status.txt
until grep -q "N11c done" $SPU/mx/night_status.txt 2>/dev/null; do sleep 30; done
cd $SPU/mx/run
timeout 1800 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); t0=tic; R1 = test_mixed_domain_pde({'convert','ops','exec'},{'C1_heat2d_scalar_w'}); R2 = test_mixed_domain_pde('exact'); fprintf('ELAPSED %.0f s\n',toc(t0)); save('$SP/N13a.mat','R1','R2');" > $SPU/mx/N13a.log 2>&1
echo "N13a done $(date)" >> $SPU/mx/night_status.txt
timeout 1800 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'D6_heat1d_reaction',struct('a0',0.9),'stab','P',true; 'S18_lib43_interior_delay',struct('tau',0.5,'sgn',-1,'a0',0.9,'a1',0),'stab','P',true; 'S18_lib43_interior_delay',struct('tau',0.5,'sgn',-1,'a0',0.9,'a1',1),'stab','P',true; 'S18_lib43_interior_delay',struct('tau',0.5,'sgn',-1,'a0',1.2,'a1',1),'stab','P',true; 'S18_lib43_interior_delay',struct('tau',0.5,'sgn',-1,'a0',1.5,'a1',1),'stab','P',true}; R = mx_exec3(c,'light'); save('$SP/N13b.mat','R');" > $SPU/mx/N13b.log 2>&1
echo "N13b done $(date)" >> $SPU/mx/night_status.txt
timeout 2400 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); [gi,h2] = mx_c1_exact(); c = {'C1_heat2d_scalar_w',struct(),'hinf','Q',gi; 'C1_heat2d_scalar_w',struct(),'h2','o',h2; 'C1_heat2d_scalar_w',struct(),'h2','oco',h2}; R = mx_exec2(c,'heavy'); save('$SP/N13c.mat','R');" > $SPU/mx/N13c.log 2>&1
echo "N13c done $(date)" >> $SPU/mx/night_status.txt
