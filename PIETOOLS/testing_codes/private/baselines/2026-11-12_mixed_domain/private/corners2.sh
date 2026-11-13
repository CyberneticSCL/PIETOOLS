#!/bin/bash
# After the 3-D run (N9): S3 and S8 at 0.9 of the margin with the corrected
# corner-row set (export t3), light.
until grep -q "N9 done" /c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx/night_status.txt 2>/dev/null; do sleep 30; done
cd /c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx/run
timeout 3600 matlab -sd "C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx/run" -batch "addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx'); mx_frozen_path('t3'); c = {'S3_edge1d_heat2d',struct('r',0.9*pi^2),'stab','P',true; 'S8_two_edges',struct('r3',0.9*pi^2),'stab','P',true}; R = mx_exec(c,'light'); save('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx/N10.mat','R');" > /c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx/N10.log 2>&1
echo "N10 done $(date)" >> /c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx/night_status.txt
