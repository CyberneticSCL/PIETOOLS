% bl_h2ctrl.m -- the control the 2-D gamma ladder needs.
%
% cuADMM called every gamma down to 0.18 feasible, 1.05x the closed-form gain
% 0.1711. That would mean bisection recovers a bound 7x tighter than the
% objective form's 1.272. Before believing it: the Hinf relaxation CANNOT be
% feasible below the true gain, so gammas underneath 0.1711 must come back
% infeasible. If they do not, the pinf<1e-3 rule is producing false positives
% and the whole ladder is an artefact of the budget rather than a result.
%
% 0.12 and 0.05 are 0.70x and 0.29x the analytic value -- unambiguously below.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
DMP  = fullfile(cuadmm_outdir(),'baseline','dumps');
BIS  = fullfile(cuadmm_outdir(),'baseline','bisect');
for g = [0.12 0.05]
    lab = sprintf('hinf2_rd_b%03.0f',g*100);
    D = bl_fix(fullfile(DMP,'hinf2_rd.mat'),g,fullfile(BIS,lab));
    fprintf('HC %.4f|%.2f|%d|%s\n',g,g/0.1711,D.m,lab);
end
fprintf('HCDONE\n');
