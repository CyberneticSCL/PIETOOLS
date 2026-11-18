% PF_RUN_C4B  Phase 3, reach: bisection of the reaction on S3 (r) and S8
% (r3) over [0.9, 1.0] of pi^2 with the N-D default storage, dw 1 pinned,
% for the scalar-rule faces, the Putinar box set and the per-direction
% faces: does the box set lose anything at the margin? Six verdicts per
% set and system. SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'facesS','box','faces'},2);
fr = [0.95 0.975 0.9875 0.99375 0.996875 0.9984375];
pf_exec('S3_edge1d_heat2d','r',fr*pi^2,C,1,{'def'},'light','c4b_s3_reach');
pf_exec('S8_two_edges','r3',fr*pi^2,C,1,{'def'},'light','c4b_s8_reach');
disp('## C4B DONE');
