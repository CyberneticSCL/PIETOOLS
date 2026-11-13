function mx_frozen_path(tree)
% MX_FROZEN_PATH(TREE) sets the MATLAB path for a run on a frozen export:
% the default path, the export TREE, MOSEK 11, the test folder of the live
% checkout (test files only) and this folder. Errors unless every probed
% name resolves once, inside the export; prints the commit line of TREE.
% INPUT
% - tree:   (optional) char, 't' (the export of 6c6f60dd, default), 't2'
%           (t with the corner-row files of set A) or 't3' (set B).
% OUTPUT
% none; the path, and appdata 'mx_tree' = TREE, on which the runners turn
% the corner rows on for 't2' and 't3'.

% Initial coding MMP, 11/12/2026: path guard of the mixed-domain runs
% (README.md); the export root is a session scratch folder.
if nargin<1, tree = 't'; end
root = ['C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/',tree,'/PIETOOLS'];
restoredefaultpath;
addpath(genpath(root));
addpath('C:/Program Files/Mosek/11.0/toolbox/r2019b');
addpath('C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/testing_codes/private/converters/tests/mixed_domain');
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/3eb3e98d-022c-4de8-a048-8544e64303d8/scratchpad/mx');
names = {'copvar','sopvar','cdopvar','PIETOOLS_stability_sop','PIETOOLS_Hinf_gain_sop', ...
         'poscopvar_direct','poscopvar','lpi_ineq_sop','get_lift_degs','convert_PIETOOLS_PDE','solve_PIETOOLS_BVP', ...
         'PIESIM','pie2piesim','apply_sopvar','lpisettings_sop','sossolve','pde_var','polynomial'};
for k=1:numel(names)
    w = which(names{k},'-all');
    w = strrep(w,'\','/');
    inroot = startsWith(w,root);
    if nnz(inroot)~=1 || ~inroot(1)
        error('mx_frozen_path:resolve','%s resolves to: %s',names{k},strjoin(w,' | '));
    end
end
assert(~isempty(which('mosekopt')),'MOSEK not on the path');
for f={'test_mixed_domain_pde','mixed_domain_systems','mixed_domain_colloc','mx_exec','mx_sys'}
    assert(numel(which(f{1},'-all'))==1,'%s not unique',f{1});
end
fprintf('frozen path OK: %s (%s)\n',root,strtrim(fileread(fullfile(root,'..','COMMIT'))));
setappdata(0,'mx_tree',tree);               % mx_exec turns the corner rows on for 't2'
end
