% PF_RUN_C2E  The half step on a symmetric box: the 2-D heat benchmark with
% Dirichlet on every edge (DD x DD, kappa* = 2 pi^2), Q at (2,1), box o1,
% box0, faces0 and facesS at rtol 1e-4, to see whether the advantage of the
% offset-0 faces over the offset-0 box factors on DD x DN (12.3360 against
% 12.3284) is tied to the Neumann side. SCRATCH, MMP, 11/15/2026.
bo = struct('rtol',1e-4,'maxsolve',24,'tmax',1800);
C = pf_cfg({'box','box0','faces0','facesS'},2);
pf_heat(2,{'DD','DD'},0,0.1,C,[1 1],1,'c2e_dddd_wq1',bo);
disp('## C2E DONE');
