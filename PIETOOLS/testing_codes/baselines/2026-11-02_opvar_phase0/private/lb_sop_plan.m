function lb_sop_plan(phase)
% LB_SOP_PLAN(PHASE), 10/28/2026: the sopvar baseline in priority order,
% into results_sop.csv (lb_full_sop).
%
% INPUT
%   phase:  'A' - light, then heavy, on every runnable row except the eight
%                 2-D rows that took over 60 s per solve on the stock path
%                 (PDE 36, 38, 39, 45, 46, 47, 50, 50/H2_norm_dual: 156 to
%                 2243 s on stock extreme)
%           'B' - veryheavy, then stripped, on the same rows; then the
%                 eight 2-D rows on light, stripped, heavy and veryheavy
%
% NOTES
%   - Resumable: rerun the same call; finished rows are skipped.

here = fileparts(mfilename('fullpath'));
S = load(fullfile(here,'catalog.mat'));  C = S.C;
C = C([C.run] & ~strcmp({C.exec},'ERROR'));
big = {'PDE:36','PDE:38','PDE:39','PDE:45','PDE:46','PDE:47','PDE:50','PDE:50/H2_norm_dual'};
rest = setdiff({C.key}, big, 'stable');
switch phase
    case 'A'
        lb_full_sop({'light','heavy'}, rest);
    case 'B'
        % No row starts after 08:00 on 10/29; rerun to resume.
        stop = '2026-10-29 08:00';
        lb_full_sop({'veryheavy','stripped'}, rest, stop);
        % The eight 2-D rows one at a time, in increasing stock extreme time
        % (156 s for PDE 47 to 2243 s for PDE 39), each on every tier
        for r = {'PDE:47','PDE:46','PDE:45','PDE:38','PDE:50/H2_norm_dual','PDE:50','PDE:36','PDE:39'}
            lb_full_sop({'light','stripped','heavy','veryheavy'}, r, stop);
        end
    otherwise
        error('lb_sop_plan:phase','phase must be ''A'' or ''B''.');
end
end
