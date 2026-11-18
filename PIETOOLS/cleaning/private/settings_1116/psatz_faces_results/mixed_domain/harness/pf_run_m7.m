% PF_RUN_M7  The two-sided bound ep I <= P <= I (PF_BOX) on S17 (r2 =
% 0.5 r2*), D4 (r2 = 0.9 pi^2) and S3 (r = 0.9 pi^2, control) at light, dw 1,
% facesS, N-D default storage, ep = 0.5, 0.1, 0.01. Prediction
% (S17_obstruction.md): S17 and D4 infeasible with a Farkas certificate at
% every ep; S3 feasible once ep is below the inverse condition number of a
% certificate. SCRATCH, MMP, 11/16/2026.
C = pf_cfg({'facesS'},2);
epl = [0.5 0.1 0.01];
pf_box('S17_transport_edge_diffusion',struct('r2',0.5*pi^2),C,1,{'def'},'light','m7_s17',epl);
pf_box('D4',0.9*pi^2,C,1,{'def'},'light','m7_d4',epl);
pf_box('S3_edge1d_heat2d',struct('r',0.9*pi^2),C,1,{'def'},'light','m7_s3',epl);
disp('## M7 DONE');
