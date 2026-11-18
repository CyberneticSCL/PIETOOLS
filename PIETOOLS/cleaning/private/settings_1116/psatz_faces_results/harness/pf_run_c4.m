% PF_RUN_C4  Phase 3: the executive path on the edge-coupled cascades S3
% (parameter r, r* = pi^2) and S8 (r3, pi^2) at 0.9 of the margin, light
% tier, every term set, dw 0 and 1, stock settings_2d storage and the N-D
% default storage (route shadow). SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'none','box','product','faces','faces0','facesS','corners'},2);
pf_exec('S3_edge1d_heat2d','r',0.9*pi^2,C,[0 1],{'stock','def'},'light','c4_s3');
pf_exec('S8_two_edges','r3',0.9*pi^2,C,[0 1],{'stock','def'},'light','c4_s8');
disp('## C4 DONE');
