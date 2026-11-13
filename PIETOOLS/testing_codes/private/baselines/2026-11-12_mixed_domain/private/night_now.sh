#!/bin/bash
# 04:30 batch: every step its own MATLAB process, serial, on frozen exports.
SPU=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad
SP=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx
LIVE=/c/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS
echo "start $(date)" > $SPU/mx/night_status.txt
cd $SPU/mx/run
# 1. final verification of the fast layers on the export of 6c6f60dd
timeout 1800 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); t0=tic; R = test_mixed_domain_pde({'convert','ops','exact','piesim'}); fprintf('ELAPSED %.0f s\n',toc(t0)); save('$SP/N1.mat','R');" > $SPU/mx/N1.log 2>&1
echo "N1 done $(date)" >> $SPU/mx/night_status.txt
# 2. the executive layer
timeout 7200 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); t0=tic; R = test_mixed_domain_pde('exec'); fprintf('ELAPSED %.0f s\n',toc(t0)); save('$SP/N2.mat','R');" > $SPU/mx/N2.log 2>&1
echo "N2 done $(date)" >> $SPU/mx/night_status.txt
# 3. reach of the edge cascades at 0.5 of the margin, light
timeout 2400 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'S3_edge1d_heat2d',struct('r',0.5*pi^2),'stab','P',true; 'S8_two_edges',struct('r3',0.5*pi^2),'stab','P',true; 'S17_transport_edge_diffusion',struct('r2',0.5*pi^2),'stab','P',true}; R = mx_exec(c,'light'); save('$SP/N3.mat','R');" > $SPU/mx/N3.log 2>&1
echo "N3 done $(date)" >> $SPU/mx/night_status.txt
# 4. corner-row variant: export + three files of the live checkout, light
rm -rf $SPU/t2 && cp -r $SPU/t $SPU/t2
for f in lpi_programming_sopvar/misc/poscopvar_direct.m lpi_programming_sopvar/misc/get_lift_degs.m lpi_programming_sopvar/lpi_ineq_sop.m; do
  cp "$LIVE/$f" "$SPU/t2/PIETOOLS/$f"
done
( echo -n "6c6f60dd + live corner-row files copied $(date +%H:%M): "; cd $SPU/t2/PIETOOLS && sha1sum lpi_programming_sopvar/misc/poscopvar_direct.m lpi_programming_sopvar/misc/get_lift_degs.m lpi_programming_sopvar/lpi_ineq_sop.m | awk '{printf "%s=%s ", $2, substr($1,1,10)}' ) > $SPU/t2/COMMIT
timeout 3600 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t2'); c = {'S3_edge1d_heat2d',struct('r',0.5*pi^2),'stab','P',true; 'S3_edge1d_heat2d',struct('r',0.9*pi^2),'stab','P',true; 'S8_two_edges',struct('r3',0.5*pi^2),'stab','P',true; 'S8_two_edges',struct('r3',0.9*pi^2),'stab','P',true; 'S17_transport_edge_diffusion',struct('r2',0.5*pi^2),'stab','P',true; 'S17_transport_edge_diffusion',struct('r2',0.9*pi^2),'stab','P',true}; R = mx_exec(c,'light'); save('$SP/N4.mat','R');" > $SPU/mx/N4.log 2>&1
echo "N4 done $(date)" >> $SPU/mx/night_status.txt
# 5. heavier tiers on the conservative cases, export of 6c6f60dd
timeout 3600 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'S3_edge1d_heat2d',struct('r',0.9*pi^2),'stab','P',true; 'S18_lib43_interior_delay',struct('tau',0.5,'sgn',-1),'stab','P',true}; R = mx_exec(c,'heavy'); save('$SP/N5.mat','R');" > $SPU/mx/N5.log 2>&1
echo "N5 done $(date)" >> $SPU/mx/night_status.txt
timeout 3600 matlab -sd "$SP/run" -batch "addpath('$SP'); mx_frozen_path('t'); c = {'S3_edge1d_heat2d',struct('r',0.9*pi^2),'stab','P',true}; R = mx_exec(c,'veryheavy'); save('$SP/N6.mat','R');" > $SPU/mx/N6.log 2>&1
echo "N6 done $(date)" >> $SPU/mx/night_status.txt
echo "all done $(date)" >> $SPU/mx/night_status.txt
