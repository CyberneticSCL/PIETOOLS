% PF_RUN_C2  Phase 2: the 2-D heat benchmark DD x DN (kappa* 12.337006),
% d = 0, ep = 0.1, R at (D,w) = (1,1), Q at (2,wQ), wQ = 1 and 2, every
% term set, certified bisection on kappa. SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'none','box','product','faces','faces0','facesS','corners'},2);
pf_heat(2,{'DD','DN'},0,0.1,C,[1 1],[1 2],'c2_ddn_d0');
disp('## C2 DONE');
