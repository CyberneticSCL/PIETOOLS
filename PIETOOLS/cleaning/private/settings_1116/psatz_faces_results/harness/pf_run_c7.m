% PF_RUN_C7  The PDE 48 estimator study: every term set at dw 0, 1, 2
% pinned, stock storage first, then the N-D default storage at dw 0 and 1.
% SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'facesS','box','box0','product','faces','faces0','corners'},2);
pf_est48(C,[0 1 2],{'stock'},'light','c7_stock');
pf_est48(C,[0 1],{'def'},'light','c7_def');
disp('## C7 DONE');
