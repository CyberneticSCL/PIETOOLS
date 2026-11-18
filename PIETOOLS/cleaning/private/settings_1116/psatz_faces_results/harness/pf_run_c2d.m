% PF_RUN_C2D  Phase 2, the half-step question: the box set and the 2^N set
% at offset 0 with Q at (2,1) (degree 2w + 2 in the acting direction,
% plain block at w), against faces0 at (2,1) (2w + 1; lo 12.3359773, hi
% 12.3367312 Farkas) and every set at (2,2) (12.3367312). rtol 1e-4.
% SCRATCH, MMP, 11/15/2026.
bo = struct('rtol',1e-4,'maxsolve',24,'tmax',1800);
C = pf_cfg({'box0','product0'},2);
pf_heat(2,{'DD','DN'},0,0.1,C,[1 1],1,'c2d_wq1_o0',bo);
disp('## C2D DONE');
