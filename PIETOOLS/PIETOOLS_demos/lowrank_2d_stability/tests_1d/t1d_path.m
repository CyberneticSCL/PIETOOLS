function D = t1d_path()                                                     % CC, 09/22/2026
% t1d_path -- path setup for the 1-D regression suite.
%
% WHAT.  Puts PIETOOLS and this directory on the path, removes the entries that
% are known to shadow it, and returns the directory names the suite needs.
% Everything is located RELATIVE TO THIS FILE, so the suite runs against the
% tree it was checked out in and cannot be aimed at another one by a stale
% working directory.
%
% WHY THIS FILE EXISTS (defect B6).  Saved pathdefs and sibling checkouts have
% repeatedly made this campaign measure a copy of the code it was not editing.
% THREE mechanisms were found live on this box while the suite was being
% written, and all three are handled here or asserted by T0:
%
%   1. NESTED WORKTREES.  pietools_path_update is a single addpath(genpath(..)),
%      and .claude/worktrees/<name>/PIETOOLS is INSIDE the tree, so genpath put
%      90 entries of a second checkout on the path -- including a copy of
%      private/bm_lm2.m that was MEASURED to differ from the current one.  It
%      happened to be ordered behind the real one, but genpath order is not an
%      invariant to rely on.  Those entries are removed below.
%   2. OUT-OF-TREE COPIES.  The scratchpad harnesses this suite descends from
%      carry same-named ancestors of almost every file here.  Any path entry
%      outside this tree that provides a name the suite owns is removed.
%   3. THE CURRENT FOLDER, which is searched BEFORE the path and cannot be
%      removed.  run_regression_1d therefore cd's here for the duration of the
%      run; T0 verifies the resolution afterwards either way.  This is not
%      hypothetical: a probe run from the scratchpad silently called the
%      scratchpad's bm_setup instead of this directory's.
%
% Removals are RETURNED and PRINTED by the caller, never silent: a shadow that
% had to be removed is information about the environment, and the point of B6
% is that it was invisible the last three times.
D.tests = fileparts(mfilename('fullpath'));
D.pkg   = fileparts(D.tests);                   % lowrank_2d_stability
D.priv  = fullfile(D.pkg,'private');            % the shipped core under test
D.root  = fileparts(fileparts(D.pkg));          % the PIETOOLS tree
pu = fullfile(D.root,'pietools_path_update.m');
if ~exist(pu,'file')
    error('t1d_path: pietools_path_update.m not found at %s',pu);
end
run(pu);
addpath(D.tests);

own = t1d_ownnames();
pp  = strsplit(path,pathsep);
kill = {};
for k = 1:numel(pp)
    p = pp{k};
    if isempty(p), continue; end
    if isKclaude(p)
        kill{end+1} = p;  continue                                  %#ok<AGROW>
    end
    if startsWith(lower([p filesep]),lower([D.root filesep])), continue; end
    % outside the tree: remove only if it provides a name the suite owns, so
    % the user's own unrelated toolboxes are left alone
    for j = 1:numel(own)
        if exist(fullfile(p,[own{j} '.m']),'file')
            kill{end+1} = p;  break                                 %#ok<AGROW>
        end
    end
end
kill = unique(kill);
if ~isempty(kill), rmpath(kill{:}); end
D.removed = kill;
D.on_path_tests = any(strcmpi(strsplit(path,pathsep),D.tests));
% One MATLAB, four threads: the house rule for this box, and it keeps the
% timings quoted in the README reproducible.
maxNumCompThreads(4);
end

function tf = isKclaude(p)
% a .claude directory anywhere in the path, i.e. a Claude Code worktree or
% artifact folder -- never the tree under test, whatever it contains
tf = any(strcmpi(strsplit(p,filesep),'.claude'));
end
