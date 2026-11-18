% PF_RUN_C6B  The library stability rows at the 'stripped' tier (the w = 1
% regime on the executive path, one raise allowed), stock storage: faces0
% (what 'auto' resolves to at w = 1), facesS, box o1, box0 and the
% per-direction faces. SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'faces0','facesS','box','box0','faces'},2);
pf_lib({'33','34','35','37','39','40'},C,{'stock'},'stripped','c6b_stripped');
disp('## C6B DONE');
