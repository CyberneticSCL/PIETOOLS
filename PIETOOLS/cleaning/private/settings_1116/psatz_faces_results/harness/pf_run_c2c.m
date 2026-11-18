% PF_RUN_C2C  Phase 2, finer resolution: the configurations tied at
% rtol 1e-3 (lo = 12.3250171 for every set at wQ = 2, and faces0 at wQ = 1)
% re-bisected at rtol 1e-4 with 24 solves and 1800 s each, to separate the
% sets at the 1e-4 level where the 10/08 record separated them (faces at
% w-1 12.336731 against the tensor set 12.326178). SCRATCH, MMP, 11/15/2026.
bo = struct('rtol',1e-4,'maxsolve',24,'tmax',1800);
C = pf_cfg({'box','facesS','faces','product','corners'},2);
pf_heat(2,{'DD','DN'},0,0.1,C,[1 1],2,'c2c_wq2_fine',bo);
C0 = pf_cfg({'faces0'},2);
pf_heat(2,{'DD','DN'},0,0.1,C0,[1 1],1,'c2c_wq1_fine',bo);
disp('## C2C DONE');
