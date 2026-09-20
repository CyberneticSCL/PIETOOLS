function pielr_path()
% PIELR_PATH  Put PIETOOLS on the MATLAB path and verify it resolves cleanly.
%
% Run this once per MATLAB session before anything else in this package.
%
% What it does, and why each step exists:
%   1. Locates PIETOOLS: if pietools_path_update is already on the path it is
%      re-run (idempotent); otherwise PIETOOLS_ROOT below must point at the
%      folder that contains pietools_path_update.m.
%   2. Strips any path entry under a '.claude' folder.  Agent worktrees under
%      such folders once put 88 shadow copies of PIETOOLS on the path, and with
%      them present EVERY 2D PDE->PIE conversion fails with the misleading
%      error "boundary conditions appear not to be invertible".  On a machine
%      with no such entries this step is a no-op.
%   3. Adds this package's folder to the path.
%   4. Asserts that every routine this package leans on resolves EXACTLY ONCE,
%      inside PIETOOLS (or inside this package for its own files), and errors
%      loudly otherwise.  A silently shadowed converter does not fail cleanly;
%      it produces wrong or misleading results, so shadowing is a hard error.
%
% ---------------------------------------------------------------------------
% EDIT THIS if PIETOOLS is not already on your MATLAB path:
PIETOOLS_ROOT = '';   % e.g. 'C:/Users/me/Documents/GitHub/PIETOOLS/PIETOOLS'
% ---------------------------------------------------------------------------

pu = which('pietools_path_update');
if isempty(pu)
    if isempty(PIETOOLS_ROOT) || ...
            ~exist(fullfile(PIETOOLS_ROOT,'pietools_path_update.m'),'file')
        error('pielr_path:noPIETOOLS', ...
            ['PIETOOLS was not found.  Either add PIETOOLS to your MATLAB path\n' ...
             '(run its pietools_path_update.m once), or open pielr_path.m and set\n' ...
             'PIETOOLS_ROOT to the folder that contains pietools_path_update.m.']);
    end
    pu = fullfile(PIETOOLS_ROOT,'pietools_path_update.m');
end
% run() cd's into the script folder, so the script's addpath(genpath(pwd))
% genpaths the PIETOOLS root regardless of the caller's cwd
run(pu);

% strip agent-worktree entries (see header step 2); harmless when none exist
p = strsplit(path,pathsep);
drop = contains(p,[filesep '.claude' filesep]);
if any(drop)
    path(strjoin(p(~drop),pathsep));
    fprintf('pielr_path: stripped %d ''.claude'' path entries\n',sum(drop));
end

here = fileparts(mfilename('fullpath'));
addpath(here);

% SeDuMi is required (certification solves a small SDP) but is NOT part of
% PIETOOLS; fail here with a clear message rather than deep inside a solve
if isempty(which('sedumi'))
    error('pielr_path:noSeDuMi', ...
          ['SeDuMi is not on the MATLAB path.  Install SeDuMi and addpath ' ...
           'it before using this package.']);
end

% single-resolution asserts.  Entries under an @class folder are method
% tables, not shadows (they only dispatch on their own class), so they are
% excluded from the count.
% free functions only: class methods (clean_opvar, convert, ...) live under
% @class folders, which the '@' filter below excludes because a method cannot
% shadow calls on other types
need = {'poslpivar_2d','lpi_eq_2d','pde_var','get_eq_opts_2D','lpiprogram', ...
        'lpigetsol','monomials','sortrows_integerTable'};
own  = {'pielr_certify','pielr_certify_pos','pielr_tensor', ...        % CC, 09/19/2026
        'set2d_deg','nb_rd2d','build_stab_2d_st2'};                    % CC, 09/19/2026
% (set2d_deg/nb_rd2d/build_stab_2d_st2 moved out of private/ because the       % CC, 09/19/2026
%  demos are now SCRIPTS, and scripts cannot see private/)                     % CC, 09/19/2026
bad = {};
for i = 1:numel(need)
    w = which(need{i},'-all');
    w = w(~cellfun(@(x)contains(x,[filesep '@']),w));
    if numel(w)~=1
        bad{end+1} = sprintf('%s (%d copies)',need{i},numel(w)); %#ok<AGROW>
    end
end
for i = 1:numel(own)
    w = which(own{i},'-all');
    w = w(~cellfun(@(x)contains(x,[filesep '@']),w));
    if numel(w)~=1 || ~strcmp(fileparts(w{1}),here)
        bad{end+1} = sprintf('%s (%d copies, expected only %s)', ...
                             own{i},numel(w),here); %#ok<AGROW>
    end
end
if ~isempty(bad)
    error(['pielr_path: shadowed or missing functions -- fix the path before ' ...
           'running anything:\n  %s'],strjoin(bad,sprintf('\n  ')));
end
fprintf('pielr_path: OK  (PIETOOLS at %s)\n',fileparts(pu));
end
