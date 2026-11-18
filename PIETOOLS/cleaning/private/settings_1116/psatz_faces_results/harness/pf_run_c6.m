% PF_RUN_C6  Regression: the 2-D library rows (stability 33, 34, 35, 37,
% 39, 40; gain 45, 46, 47; estimator 48) at light with the tier's own raise
% loop, stock storage, four term sets: the scalar-rule faces (today's
% 'auto'), the Putinar box set, the per-direction faces and the 2^N set.
% SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'facesS','box','faces','product'},2);
pf_lib({'33','34','35','37','39','40','45','46','47','48'},C,{'stock'},'light','c6_lib2d');
disp('## C6 DONE');
