function info = p0_path()
% INFO = P0_PATH() sets the MATLAB path for the Phase-0 fixtures: the frozen
% export of the recorded commit, MOSEK and SeDuMi, and nothing of the main
% checkout, which other sessions edit while the fixtures are recorded.
%
% OUTPUT
% - info:   'struct' with fields
%             root     char, the export the path points to;
%             commit   char, the commit the export was made from;
%             git      char, the same, in the form lb_run logs;
%             version  char, MATLAB version;
%             mosek    char, MOSEK version, or 'none';
%             threads  double, maxNumCompThreads.
%
% NOTES
% The export is 'git archive <commit>:PIETOOLS'; its root holds
% P0_COMMIT.txt. The environment variable P0_ROOT names the export (passes
% 1-2: p0_eb5d0600, passes 3-4: p0_c81c3d04, the default). An error is
% raised when a name used by the fixtures resolves to more than one file,
% or first to a file outside the export.
ROOT = getenv('P0_ROOT');
if isempty(ROOT),   ROOT = 'C:\Users\mpeet\wt\p0_c81c3d04';     end
if ~exist(fullfile(ROOT,'P0_COMMIT.txt'),'file')
    error('p0_path:root','No export at %s (README, Tree, recreates it).',ROOT);
end
restoredefaultpath;
addpath(genpath(ROOT));
warning('off','MATLAB:rmpath:DirNotFound');
% Research harnesses with common function names; never used by the fixtures
for f = {'PIETOOLS_demos\cuadmm','PIETOOLS_demos\lowrank_2d_stability','cleaning'}
    rmpath(genpath(fullfile(ROOT,f{1})));
end
addpath('C:\Program Files\Mosek\11.0\toolbox\r2019b');
sedumi = 'C:\YALMIP\yalmip_install\tbxmanager\toolboxes\sedumi';
if exist(sedumi,'dir'),     addpath(genpath(sedumi));   end
names = {'convert','convert_PIETOOLS_PDE','convert_PIETOOLS_PDE_sopvar','PIESIM', ...
         'PIESIM_time_integrate','pie2opvar','op2copvar_sop','apply_sopvar','apply_opvar', ...
         'lpiprogram','lpidecvar','lpisetobj','lpisolve','sossolve','lpisettings_sop', ...
         'PIETOOLS_stability_sop','PIETOOLS_Hinf_gain_sop','lpi_ineq_sop','poscopvar', ...
         'examples_pde_library_PIESIM_1D','examples_pde_library_PIESIM_2D', ...
         'DEMO1_stability_and_simulation'};
for k = 1:numel(names)
    w = which(names{k},'-all');
    w = w(endsWith(w,'.m') & ~contains(w,'Shadowed'));
    if strcmp(names{k},'convert')
        ok = any(contains(w,'@pde_struct')) && all(startsWith(w,ROOT));
    else
        ok = numel(w)==1 && startsWith(w{1},ROOT);
    end
    if ~ok
        error('p0_path:resolve','%s resolves to %d files:\n%s',names{k},numel(w),strjoin(w,newline));
    end
end
info.root = ROOT;
info.commit = strtrim(fileread(fullfile(ROOT,'P0_COMMIT.txt')));
info.git = [info.commit(1:8),' (export, clean)'];
info.version = version;
try
    [~,res] = mosekopt('version');
    info.mosek = sprintf('%d.%d.%d',res.version.major,res.version.minor,res.version.revision);
catch
    info.mosek = 'none';
end
info.threads = maxNumCompThreads;
end
