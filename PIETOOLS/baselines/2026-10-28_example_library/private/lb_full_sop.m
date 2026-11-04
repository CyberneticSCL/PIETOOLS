function lb_full_sop(tiers, ids, stopat)
% LB_FULL_SOP(TIERS, IDS) runs the library baseline on the sopvar path:
% each catalog row through the _sop executive and form that
% PIETOOLS_auto_execute_sop uses for its stock executive, with
% lpisettings_sop(tier,0,'','mosek'), into results_sop.csv, one tier at a
% time in the order given.
%
% INPUT
%   tiers:  optional cellstr of settings_sopvar tiers (default {'light',
%           'heavy','veryheavy','stripped'}: the default tier first;
%           custom is light in settings_sopvar and is not run)
%   ids:    optional catalog selection as in lb_run (default 'all')
%   stopat: optional date-time string; no row starts after it (lb_run)
%
% NOTES
%   - Resumable as lb_run: rerun the same call after an interruption.
%   - PAUSE_sop in this folder stops the run before the next row.

if nargin<1 || isempty(tiers), tiers = {'light','heavy','veryheavy','stripped'}; end
if nargin<2 || isempty(ids), ids = 'all'; end
here = fileparts(mfilename('fullpath'));
if nargin<3, stopat = ''; end
opts = struct('path','sop','results',fullfile(here,'results_sop.csv'),'max2d',600,'stopat',stopat);
for t = 1:numel(tiers)
    if exist(fullfile(here,'PAUSE_sop'),'file'), return, end
    lb_run(ids, tiers(t), opts);
end
end
