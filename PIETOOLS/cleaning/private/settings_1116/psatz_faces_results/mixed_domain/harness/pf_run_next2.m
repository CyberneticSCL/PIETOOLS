% PF_RUN_NEXT2  M7 (two-sided bounds on P: S17, D4, S3) then M8 (S18
% attribution and soundness), one MATLAB process. SCRATCH, MMP, 11/16/2026.
try,    run('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces/pf_run_m7.m');
catch ME,   fprintf('@@ ERROR m7: %s\n',strrep(ME.message,newline,' '));  end
try,    run('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces/pf_run_m8.m');
catch ME,   fprintf('@@ ERROR m8: %s\n',strrep(ME.message,newline,' '));  end
disp('## NEXT2 DONE');
