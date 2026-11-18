% PF_RUN_C1B  Phase 1b: the converter-built pinned-P plants of E2 (local
% coefficients L1, L2 and the ODE-coupled plant), absolute lam rungs as in
% the 09/24 record, every term set at dw 0 and 1. SCRATCH, MMP, 11/15/2026.
rungs = [0.10 0.50 0.90 0.99 0.999 0.9999];
pf_pinned(2,{'DD','DD'},{'prodonly','prodonly0'},[0 1 2],rungs,'c1_prodonly');   % the stock product alone
sets = {'none','box','product','faces','facesS','corners'};
lams = [1 2 5 10 15 19];
pf_pinned_pde('L1',sets,[0 1],lams,'l1');
pf_pinned_pde('L2',sets,[0 1],lams,'l2');
pf_pinned_pde('ODE',sets,[0 1],lams,'ode');
disp('## C1B DONE');
