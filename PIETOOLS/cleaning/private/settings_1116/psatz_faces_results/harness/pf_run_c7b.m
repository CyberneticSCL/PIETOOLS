% PF_RUN_C7B  PDE 48 estimator, the ladder one step up: faces0, facesS,
% box and box0 at dw 3 (stock storage), against the record's best bound
% 0.0024154 at 2.41M variables and faces0's 0.0045059 at dw 1 (110k).
% SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'faces0','facesS','box','box0'},2);
pf_est48(C,3,{'stock'},'light','c7b_dw3');
disp('## C7B DONE');
