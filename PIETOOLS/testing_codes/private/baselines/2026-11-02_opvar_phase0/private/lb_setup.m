function info = lb_setup()
% INFO = LB_SETUP() puts PIETOOLS and MOSEK on the path and asserts that the
% functions the library baseline calls resolve to exactly one file inside
% the PIETOOLS root.
%
% OUTPUT
%   info:   struct with fields 'root' (PIETOOLS root folder), 'mosek'
%           (path of mosekopt), 'version' (MATLAB version), 'threads'
%           (maxNumCompThreads), 'git' (HEAD hash and dirty files of the
%           repository, read-only), 'here' (this folder).
%
% NOTES
%   - restoredefaultpath first, so a saved pathdef cannot shadow toolbox
%     files; heatNd_path then runs pietools_path_update, strips nested
%     worktrees and asserts single resolution of the core functions.
%   - The scratch folder stays first on the path through cd(here).
%   - Phase 0, 11/02/2026: the path is the frozen export of P0_PATH, not the
%     main checkout, which other sessions edit during the runs; info.git
%     names the export's commit.
%   - MMP, 11/08/2026: the stock executives and the stock settings left the
%     list of functions asserted to resolve (BEGIN/END block below): they
%     are retired, and the 'sop' path of LB_RUN, the one P0_F3 runs, calls
%     neither; with them in the list every F3 run stopped here.

here = fileparts(mfilename('fullpath'));
% restoredefaultpath;                                                       % Phase 0, 11/02/2026 (was)
% addpath('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS\PIETOOLS_demos\sopvar_demos'); % Phase 0, 11/02/2026 (was)
% hp = heatNd_path();                                                       % Phase 0, 11/02/2026 (was)
hp = p0_path();                                 % frozen export             % Phase 0, 11/02/2026
cd(here);
info = hp;
info.here = here;

% Functions of the stock flow, not covered by heatNd_path: one file each.
% BEGIN MMP, 11/08/2026: the stock executives (retired to cleaning/private,
% f83486c8) and the stock settings lpisettings, settings_PIETOOLS_light and
% settings_PIETOOLS_light_2D (leaving the path) are no longer asserted. The
% list was:
% fns = {'convert','convert_PIETOOLS_PDE','initialize_PIETOOLS_PDE', ...
%        'convert_PIETOOLS_DDE','initialize_PIETOOLS_DDE', ...
%        'convert_PIETOOLS_NDS','initialize_PIETOOLS_NDS', ...
%        'convert_PIETOOLS_DDF','initialize_PIETOOLS_DDF','minimize_PIETOOLS_DDF', ...
%        'lpisettings','settings_PIETOOLS_light','settings_PIETOOLS_light_2D', ...
%        'PIETOOLS_PIE2PDEstability','PIETOOLS_PIE2PDEstability_dual', ...
%        'PIETOOLS_Hinf_gain','PIETOOLS_Hinf_gain_dual','PIETOOLS_Hinf_estimator', ...
%        'PIETOOLS_Hinf_control','PIETOOLS_H2_norm_o','PIETOOLS_H2_norm_c', ...
%        'PIETOOLS_H2_estimator','PIETOOLS_H2_control', ...
%        'PIETOOLS_PIE2PDEstability_2D','PIETOOLS_Hinf_gain_2D', ...
%        'PIETOOLS_PDE_Ex_Reaction_Diffusion_Gahlawat', ...
%        'lpisettings_sop','pie2copvar','PIETOOLS_stability_sop', ...
%        'PIETOOLS_Hinf_gain_sop','PIETOOLS_H2_norm_sop','PIETOOLS_Hinf_estimator_sop', ...
%        'PIETOOLS_Hinf_control_sop','PIETOOLS_H2_estimator_sop','PIETOOLS_H2_control_sop'};
fns = {'convert','convert_PIETOOLS_PDE','initialize_PIETOOLS_PDE', ...
       'convert_PIETOOLS_DDE','initialize_PIETOOLS_DDE', ...
       'convert_PIETOOLS_NDS','initialize_PIETOOLS_NDS', ...
       'convert_PIETOOLS_DDF','initialize_PIETOOLS_DDF','minimize_PIETOOLS_DDF', ...
       'PIETOOLS_PDE_Ex_Reaction_Diffusion_Gahlawat', ...
       'lpisettings_sop','pie2copvar','PIETOOLS_stability_sop', ...
       'PIETOOLS_Hinf_gain_sop','PIETOOLS_H2_norm_sop','PIETOOLS_Hinf_estimator_sop', ...
       'PIETOOLS_Hinf_control_sop','PIETOOLS_H2_estimator_sop','PIETOOLS_H2_control_sop'};
% END MMP, 11/08/2026
for k = 1:numel(fns)
    w = which(fns{k},'-all');
    w = w(endsWith(w,'.m') & ~contains(w,'Shadowed'));
    inroot = startsWith(strrep(w,'\','/'),strrep(hp.root,'\','/'));
    % 'convert' also has class methods (@pde_struct and others): count
    % only the files outside class folders, plus the pde_struct method.
    if strcmp(fns{k},'convert')
        ok = any(contains(w,'@pde_struct'));
    else
        ok = nnz(inroot & ~contains(w,'@'))==1 && inroot(1);
    end
    if ~ok
        error('lb_setup:resolve','%s resolves to %d files:\n%s', ...
              fns{k},numel(w),strjoin(w,newline));
    end
end

% Repository state, read-only: the baseline depends on the working tree.
% [~,h] = system(['git -C "',hp.root,'" rev-parse --short HEAD']);          % Phase 0, 11/02/2026 (was)
% [~,s] = system(['git -C "',hp.root,'" status --short --untracked-files=no']); % Phase 0, 11/02/2026 (was)
% info.git = strtrim([strtrim(h),' dirty: ',strjoin(strsplit(strtrim(s),newline),' | ')]); % Phase 0, 11/02/2026 (was)
info.git = hp.git;                              % the export has no .git    % Phase 0, 11/02/2026
end
