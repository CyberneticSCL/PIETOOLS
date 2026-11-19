function lb_ddeq4()
% LB_DDEQ4() the DDE and NDSDDF library rows under the light tier with
% deg_Q{1} = [4 2 2] (maintainer, 11/18/2026 01:28), into results_ddeQ4.csv
% of the 11/17 harness folder (the current folder); size gate as
% LB_TIERS_1117. MMP, 11/18/2026.
here = pwd;
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
opts = struct('path','sop','results',fullfile(here,'results_ddeQ4.csv'),'max2d',600);
lb_run('DDE', {'light'}, opts);
lb_run('NDSDDF', {'light'}, opts);
fprintf('## LB_DDEQ4 DONE %s\n',char(datetime('now')));
end
