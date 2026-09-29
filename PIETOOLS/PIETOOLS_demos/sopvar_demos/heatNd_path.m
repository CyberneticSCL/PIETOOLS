function info = heatNd_path()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = HEATND_PATH() puts this checkout of PIETOOLS on the path for the
% heatNd benchmark and asserts that every routine the benchmark relies on
% resolves exactly once, inside the checkout.
%
% Why: a saved pathdef or a nested worktree has silently run the wrong
% copy before, and same-named files ahead on the path crashed MATLAB.
% Steps: pietools_path_update (genpath of the checkout); drop entries with
% '.claude' or 'worktrees'; add Mosek's toolbox if it is not on the path;
% assert single in-checkout resolution (class methods of one name are one
% owner each, not shadows).
%
% OUTPUT struct: threads (MATLAB's maxNumCompThreads; MOSEK does NOT
% inherit it - measured 24 MOSEK threads with MATLAB at 8; HEATND_SOLVE
% records MOSEK's own count), version, mosek (path of mosekopt), root.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): corrected the thread comment.
% CC, 09/28/2026 (move): folder moved from sopvar/Testfolder/sdopvar/
%                  claude_tests/heatNd to PIETOOLS_demos/sopvar_demos, so
%                  ROOT is two levels up, not five.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
% ROOT = fileparts(fileparts(fileparts(fileparts(fileparts(here)))));   % .../PIETOOLS/PIETOOLS  % CC, 09/28/2026 (was)
ROOT = fileparts(fileparts(here));          % .../PIETOOLS/PIETOOLS         % CC, 09/28/2026
run(fullfile(ROOT,'pietools_path_update.m'));
p = strsplit(path,pathsep);
bad = contains(p,'.claude') | contains(p,'worktrees');
if any(bad),    path(strjoin(p(~bad),pathsep));     end
MOSEK = 'C:\Program Files\Mosek\11.0\toolbox\r2019b';
if isempty(which('mosekopt')) && exist(MOSEK,'dir'),    addpath(MOSEK);     end
fns = {'sopvar','sdopvar','copvar','cdopvar','lpivar_cdopvar','poscopvar', ...
       'copquadvar','lpi_eq_cdopvar','lpi_eq_sdopvar','eq_opts_sopvar', ...
       'degbalance_core','int_semisep','sosquadvar','dpvar2sdvar', ...
       'sossolve','lpisolve','sosprogram','lpidecvar','lpiprogram', ...
       'spblkdiag','spantiblkdiag','Sedumi2Mosek','MosekSol2SedumiSol', ...
       'pde_var','cx_resid'};
d = dir(fullfile(here,'heatNd_*.m'));
fns = [fns, cellfun(@(f) f(1:end-2),{d.name},'UniformOutput',false)];
for k = 1:numel(fns)
    w = which(fns{k},'-all');
    w = w(endsWith(w,'.m') & ~contains(w,'Shadowed'));
    % One hit per owner (class folder, or none for a plain function): a
    % class method of the same name is dispatch, a second plain file or a
    % second copy of a class is a shadow.
    % A MATLAB toolbox file of the same name BEHIND the checkout's copy is
    % not reached (e.g. mbc's spblkdiag); one AHEAD of it would be.
    own = regexp(w,'@(\w+)[\\/]','tokens','once');
    own = cellfun(@(z) strjoin(z,''),own,'UniformOutput',false);
    inroot = startsWith(strrep(w,'\','/'),strrep(ROOT,'\','/'));
    if isempty(w) || ~inroot(1) || numel(unique(own(inroot)))~=nnz(inroot)
        error('heatNd_path:resolve','%s resolves %d times / outside %s:\n%s', ...
              fns{k},numel(w),ROOT,strjoin(w,newline));
    end
end
if isempty(which('mosekopt')),  warning('heatNd_path:mosek','mosekopt not on path');  end
info = struct('threads',maxNumCompThreads,'version',version, ...
              'mosek',which('mosekopt'),'root',ROOT);
end
