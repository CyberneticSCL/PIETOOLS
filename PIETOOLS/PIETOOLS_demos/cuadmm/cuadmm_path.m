function cuadmm_path(mode)                                                  % CC, 09/25/2026
% CUADMM_PATH  Put THIS checkout's PIETOOLS and the SDP solvers on the path, and
% verify that everything this package depends on resolves exactly once, inside
% this checkout.
%
%   cuadmm_path            % non-destructive: adds to the current path
%   cuadmm_path('reset')   % restoredefaultpath first, as the 2026-09-24
%                          % baseline was measured
%
% Every script in this package calls it first.  It replaces the scratchpad's
% pp.m, which hardcoded this workstation's paths and always ran
% restoredefaultpath -- acceptable in a throwaway -batch process, not in
% someone's interactive session, so resetting is now opt-in.
%
% ANCHORED TO THE FILE, NOT TO THE PATH.  ROOT is the PIETOOLS folder this
% cuadmm_path.m lives in (three levels up).  Everything below is relative to it:
% which copy to put on the path, which entries to strip, and where every
% dependency must resolve.  A path-relative version had two measured-by-review
% failures: 'reset' could not find pietools_path_update once the path was wiped,
% and stripping every '.claude'/'worktrees' entry would delete the copy UNDER
% TEST when run from an agent worktree, leaving the main copy to satisfy the
% resolution check -- the "tests the wrong copy" failure the check exists for.
%
% Steps:
%   1. Removes any cuadmm_shadow folders, with a warning.  A shadow replaces
%      lpisolve or poslpivar for EVERY caller while present, so a leftover
%      from an earlier script must not survive into the next run.  Note this
%      also removes a shadow a caller installed on purpose -- so never call a
%      harness runner (bl_run, bl_check, ...) with a shadow active.
%   2. Runs ROOT/pietools_path_update.m.
%   3. Strips entries under '.claude' or 'worktrees' that are NOT inside ROOT
%      -- other checkouts, never this one.
%   4. Adds SeDuMi and Mosek if the folders below exist.
%   5. ASSERTS each routine resolves exactly once AND inside ROOT.

% ---------------------------------------------------------------------------
% EDIT THESE for your machine. Missing folders are skipped.
SEDUMI = 'C:\Users\mpeet\ASU Dropbox\Matthew Peet\Codes\setup_nonlinear_Matlab\SeDuMi_1_3';
MOSEK  = 'C:\Program Files\Mosek\11.0\toolbox\r2019b';
% ---------------------------------------------------------------------------

ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
pu   = fullfile(ROOT,'pietools_path_update.m');
if ~exist(pu,'file')
    error('cuadmm_path:noPIETOOLS', ...
        'expected pietools_path_update.m at %s -- is this package inside PIETOOLS_demos?',pu);
end

if nargin >= 1 && strcmpi(mode,'reset'), restoredefaultpath; end

% 1. stale shadows
p = strsplit(path,pathsep);
sh = contains(p,'cuadmm_shadow_');
if any(sh)
    path(strjoin(p(~sh),pathsep));
    warning('cuadmm_path:shadowRemoved','removed %d active cuadmm_shadow folder(s)',sum(sh));
end

% 2. this checkout's PIETOOLS. run() cd's into the script's folder, so its
% genpath covers ROOT regardless of the caller's working directory.
run(pu);

% 3. worktree copies. The pattern is tested on the part of each entry AFTER
% ROOT, so two cases are both handled:
%   - ROOT itself inside a worktree (running from .claude/worktrees/x/PIETOOLS):
%     its own entries have no '.claude' after ROOT, so they are KEPT;
%   - worktrees NESTED inside this checkout (ROOT/.claude/worktrees/...), which
%     pietools_path_update's genpath picks up: stripped. Measured 2026-09-25:
%     keeping everything under ROOT left ROOT/.claude/worktrees/modest-bassi-946e4a
%     on the path and poslpivar resolved twice.
% Entries outside ROOT are tested whole.
p = strsplit(path,pathsep);
mine = strcmpi(p,ROOT) | strncmpi(p,[ROOT filesep],numel(ROOT)+1);   % not a sibling like PIETOOLS_old
rest = p;
rest(mine) = cellfun(@(s) s(numel(ROOT)+1:end), p(mine), 'UniformOutput', false);
bad  = ~cellfun(@isempty, regexpi(rest,'[\\/]\.claude([\\/]|$)|[\\/]worktrees([\\/]|$)','once'));
if any(bad), path(strjoin(p(~bad),pathsep)); end

% 4. solvers
if exist(SEDUMI,'dir'), addpath(SEDUMI); end
if exist(MOSEK,'dir'),  addpath(MOSEK);  end

% 5. single resolution, inside ROOT. 'monomials' is excluded deliberately: its
% second hit is @xregcubic/monomials, a class method that cannot shadow.
need = {'poslpivar','lpi_eq','lpiprogram','lpisolve','sossolve', ...
        'poslpivar_2d','lpi_eq_2d', ...
        'bl_run','bl_check','cuadmm_settings','cuimport','bl_fix'};
for i = 1:numel(need)
    w = unique(which(need{i},'-all'));
    if numel(w) ~= 1
        error('cuadmm_path:shadow','%s resolves %d times:\n%s',need{i},numel(w), ...
              strjoin(w,newline));
    end
    if ~strncmpi(w{1},[ROOT filesep],numel(ROOT)+1)
        error('cuadmm_path:wrongCopy','%s resolves OUTSIDE this checkout:\n  %s\n  (ROOT = %s)', ...
              need{i},w{1},ROOT);
    end
end
end
