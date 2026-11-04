function p0_f3(pass, tiers)
% P0_F3(PASS,TIERS) records fixture F3: the container (_sop) executive rows
% of the 10/28 example-library catalog on the frozen export, through the
% copied 10/28 harness (lb_run with path 'sop').
%
% INPUT
% - pass:   double, the pass number; output ../F3/pass<PASS>/sop_results.csv;
% - tiers:  cell of settings names, e.g. {'light','heavy'}.
%
% NOTES
% - Columns as in the 10/28 README; the note of each row now ends with the
%   executive's classifier status, SDP shape (m, nx, ndv, Ks) and witness
%   (P0_INFONOTE through the patched lb_exec).
% - PDE 36 and 38 are skipped on both tiers: on 10/28 they took 1,653 s and
%   1,638 s at light (infeasible) and were skipped at heavy by projection;
%   with no light row the projection would run them at heavy.
% - PDE 50 and 50/H2_norm_dual are skipped on every tier since 11/02, 06:30:
%   an error on 10/28, on the export the light primal took 2,421 s (ndv
%   5.08e6, m 40,881, status unknown). The pass 1 light primal row, run
%   before the skip, stays in the results.
% - Rows already in the results file are not rerun (lb_run resumes).
here = fileparts(mfilename('fullpath'));
out = fullfile(fileparts(here),'F3',sprintf('pass%d',pass));
if ~exist(out,'dir'),   mkdir(out);     end
skip = {};
for t = {'light','heavy','veryheavy','stripped'}
    skip = [skip, {['PDE:36@',t{1}], ['PDE:38@',t{1}], ...
                   ['PDE:50@',t{1}], ['PDE:50/H2_norm_dual@',t{1}]}];        %#ok<AGROW>
end
opts = struct('path','sop','results',fullfile(out,'sop_results.csv'),'max2d',600,'skip',{skip});
for t = 1:numel(tiers)
    fprintf('\n@@@ F3 pass %d tier %s %s\n', pass, tiers{t}, char(datetime('now')));
    lb_run('all', tiers(t), opts);
end
fprintf('\nDONE F3 pass %d\n', pass);
end
