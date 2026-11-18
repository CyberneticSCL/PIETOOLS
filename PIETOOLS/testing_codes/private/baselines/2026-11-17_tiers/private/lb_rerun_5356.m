function lb_rerun_5356()
% LB_RERUN_5356() reruns PDE 53-56 (the 3-D cases and the 2-D delayed heat
% equation) under light and heavy into results_A.csv, after their error
% rows of the 11/17/2026 pass (lb_build_system before its try/catch guard)
% have been deleted from the file (lb_run resumes by key). Size gate as
% LB_TIERS_1117. MMP, 11/17/2026.
here = fileparts(mfilename('fullpath'));
cd(here);
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
opts = struct('path','sop','results',fullfile(here,'results_A.csv'),'max2d',600);
lb_run({'PDE:53','PDE:54','PDE:55','PDE:56'}, {'light','heavy'}, opts);
fprintf('## LB_RERUN_5356 DONE %s\n',char(datetime('now')));
end
