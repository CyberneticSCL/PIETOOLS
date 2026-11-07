% bl_h2probe.m -- how far down does the 2-D Hinf bound go under BISECTION?
%
% The objective form returned gam = 1.272 on hinf2_rd, against the shipped
% example's closed-form L2 gain of 0.1711 -- 7.4x conservative, with Mosek
% reporting UNKNOWN. But the pinned feasibility form is already satisfiable at
% 0.80*1.272 = 1.018 (cuADMM, pinf 2.45e-05 at 5000 iterations), so the
% objective form did not find the best bound its own relaxation admits.
%
% This writes the pinned programs on a ladder running from just above the
% analytic value up to the objective form's answer. If the relaxation is
% feasible well below 1.272, bisection recovers a materially tighter bound than
% minimising gamma did -- which is the practical case for switching the suite to
% bisection, independently of which solver runs the step.
%
% Mosek is not run here: it took 341-480 s per step on this case and returned
% UNKNOWN every time, so it cannot adjudicate the ladder. The dumps are written
% for cuADMM, whose per-step cost at m=14006 was 172 s for 5000 iterations.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
DMP  = fullfile(cuadmm_outdir(),'baseline','dumps');
BIS  = fullfile(cuadmm_outdir(),'baseline','bisect');
GEX  = 0.1711;                      % closed-form gain from the example file
GAMS = [0.18 0.25 0.40 0.60 0.85];  % between the analytic value and 1.018
fprintf('HP gam|gam_over_exact|m|wrote\n');
for g = GAMS
    lab = sprintf('hinf2_rd_a%03.0f',g*100);
    D = bl_fix(fullfile(DMP,'hinf2_rd.mat'),g,fullfile(BIS,lab));
    fprintf('HP %.4f|%.2f|%d|%s\n',g,g/GEX,D.m,lab);
end
fprintf('HPDONE\n');
