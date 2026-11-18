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
%
% MMP, 11/10/2026: heatNd_path and cx_resid, which it asserts, moved to
%                  testing_codes/private/; both folders are added in place
%                  of PIETOOLS_demos/sopvar_demos.

here = fileparts(mfilename('fullpath'));
restoredefaultpath;
% addpath('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS\PIETOOLS_demos\sopvar_demos');  % MMP, 11/10/2026 (was)
% heatNd_path and cx_resid: in testing_codes/private/ since 11/10/2026.     % MMP, 11/10/2026
% addpath('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS\testing_codes\private\PIETOOLS_demos\sopvar_demos', ...  % MMP, 11/17/2026 (was)
%         'C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS\testing_codes\private\sopvar\Testfolder\sdopvar\claude_tests\cx_exec');  % MMP, 11/10/2026 (was)
% The checkout that holds this harness, so a copy in a worktree runs that   % MMP, 11/17/2026
% worktree's tree (HEAD reference run of the tier pass, 11/17/2026)         % MMP, 11/17/2026
root = fileparts(fileparts(fileparts(fileparts(fileparts(here)))));         % MMP, 11/17/2026
addpath(fullfile(root,'testing_codes','private','PIETOOLS_demos','sopvar_demos'), ...  % MMP, 11/17/2026
        fullfile(root,'testing_codes','private','sopvar','Testfolder','sdopvar','claude_tests','cx_exec'));  % MMP, 11/17/2026
hp = heatNd_path();
cd(here);
info = hp;
info.here = here;

% Functions of the stock flow, not covered by heatNd_path: one file each.
fns = {'convert','convert_PIETOOLS_PDE','initialize_PIETOOLS_PDE', ...
       'convert_PIETOOLS_DDE','initialize_PIETOOLS_DDE', ...
       'convert_PIETOOLS_NDS','initialize_PIETOOLS_NDS', ...
       'convert_PIETOOLS_DDF','initialize_PIETOOLS_DDF','minimize_PIETOOLS_DDF', ...
       ... % MMP, 11/17/2026: the stock settings and executives (lpisettings, settings_PIETOOLS_*, PIETOOLS_PIE2PDEstability*, PIETOOLS_Hinf_*, PIETOOLS_H2_*, *_2D) left the path 11/08/2026 (cleaning/private); only the sopvar path is checked
       'PIETOOLS_PDE_Ex_Reaction_Diffusion_Gahlawat', ...
       'lpisettings_sop','pie2copvar','PIETOOLS_stability_sop', ...
       'PIETOOLS_Hinf_gain_sop','PIETOOLS_H2_norm_sop','PIETOOLS_Hinf_estimator_sop', ...
       'PIETOOLS_Hinf_control_sop','PIETOOLS_H2_estimator_sop','PIETOOLS_H2_control_sop'};
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
[~,h] = system(['git -C "',hp.root,'" rev-parse --short HEAD']);
[~,s] = system(['git -C "',hp.root,'" status --short --untracked-files=no']);
info.git = strtrim([strtrim(h),' dirty: ',strjoin(strsplit(strtrim(s),newline),' | ')]);
end
