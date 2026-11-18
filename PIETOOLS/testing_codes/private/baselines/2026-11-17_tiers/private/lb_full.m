function lb_full(plan, stopat, dry)
% LB_FULL(PLAN, STOPAT, DRY) launches the full library baseline: every
% catalog row with run=yes, on all six settings tiers, into results.csv.
%
% INPUT
%   plan:   'reduced' (default) - all rows on all tiers except DDE16
%                                 (Hinf control, PIE size [3 9]) on
%                                 veryheavy and custom, projected at
%                                 2.0 h and 1.4 h per solve
%           'full'              - every row on every tier
%   stopat: optional date-time string; no row starts after it, e.g.
%           '2026-10-29 10:30' (default: no limit)
%   dry:    optional, true lists the pending rows and solves nothing
%
% NOTES
%   - Resumable: rerun the same call after an interruption or a MATLAB
%     crash; finished (library, id, settings) rows are skipped.
%   - 2-D rows: lb_run skips a row whose projected time per solve exceeds
%     600 s, projecting from its lighter tiers with a growth of at least 1 per
%     step into stripped or light, 10 into custom, heavy or veryheavy
%     (measured on PDE 35).

if nargin<1 || isempty(plan), plan = 'reduced'; end
if nargin<2, stopat = ''; end
if nargin<3, dry = false; end
opts = struct('max2d',600,'stopat',stopat,'dry',dry);
switch plan
    case 'reduced'
        opts.skip = {'DDE:DDE16@veryheavy','DDE:DDE16@custom'};
    case 'full'
        opts.skip = {};
    otherwise
        error('lb_full:plan','plan must be ''reduced'' or ''full''.');
end
lb_run('all','all',opts);
end
