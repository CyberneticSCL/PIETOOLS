function d = cuadmm_outdir()                                                % CC, 09/25/2026
% CUADMM_OUTDIR  Root folder for everything this package GENERATES: SDP dumps,
% cuADMM problem directories and logs, result tables, check output.
%
% Deliberately outside the repository by default. One overnight baseline wrote
% ~1.9 GB of dumps, which must never land in a git working tree.
%
%   setenv('CUADMM_OUT','D:\cuadmm')   % choose a location
%   d = cuadmm_outdir()                % default: <tempdir>/cuadmm
%
% Also creates baseline/ and baseline/dumps/ under it. The runners open files
% there directly; before the move they sat in a folder that already existed, so
% without this bl_check's very first fopen fails on a fresh machine.
%
% RESUME STATE LIVES HERE AND PERSISTS ACROSS SESSIONS. bl_run, bl_scale and
% bl_verify skip cases already recorded 'ok' in this folder. After changing the
% code under test, point CUADMM_OUT at a fresh folder or delete the tables, or
% the new run will silently reuse rows measured on the old code.
%
% Inputs that belong to the suite -- bl_expect.tsv, the banked expectations,
% and the frozen tables in results/ -- are read from the suite folder,
% testing_codes/private/PIETOOLS_demos/cuadmm, not from here.

% CC, 11/08/2026: help text corrected; the suite (bl_check, bl_expect.tsv,
%   results/) moved to testing_codes/PIETOOLS_demos/cuadmm (a9e6487c), then
%   to testing_codes/private/PIETOOLS_demos/cuadmm (5388c808), while this
%   file stays in PIETOOLS_demos/cuadmm with the solver bridge. Was:
%   "Inputs that belong to the package -- ... -- are read from the package,
%   not from here."
d = getenv('CUADMM_OUT');
if isempty(d), d = fullfile(tempdir,'cuadmm'); end
if ~exist(fullfile(d,'baseline','dumps'),'dir'), mkdir(fullfile(d,'baseline','dumps')); end
end
