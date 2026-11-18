% PF_RUN_NEXT4  the Q-form batch, then the storage-cone batch, one MATLAB
% process (the fixed-operator batch was cut after S18_Pstar: every target
% certified under the reader). SCRATCH, MMP, 11/16/2026.
try,    run('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces/pf_run_q.m');
catch ME,   fprintf('@@ ERROR q: %s\n',strrep(ME.message,newline,' '));  end
try,    run('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces/pf_run_cone.m');
catch ME,   fprintf('@@ ERROR cone: %s\n',strrep(ME.message,newline,' '));  end
disp('## NEXT4 DONE');
