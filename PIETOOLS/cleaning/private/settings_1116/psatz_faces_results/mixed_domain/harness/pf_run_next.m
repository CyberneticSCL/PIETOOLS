% PF_RUN_NEXT  chk2 (both certificate families) then M6 (S18 with the delay
% direction split at raised degrees), one MATLAB process. SCRATCH, MMP, 11/16/2026.
try,    run('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces/pf_run_chk2.m');
catch ME,   fprintf('@@ ERROR chk2: %s\n',strrep(ME.message,newline,' '));  end
try,    run('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces/pf_run_m6.m');
catch ME,   fprintf('@@ ERROR m6: %s\n',strrep(ME.message,newline,' '));  end
disp('## NEXT DONE');
