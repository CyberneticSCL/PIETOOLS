% PF_RUN_C01  Phase 1: the pinned-P rung sweeps, C0 (1-D DD) and C1 (2-D
% DD x DD), every set at dw 0 and 1 (dw 2 in a later job where a set has
% not reached the top rung). SCRATCH, MMP, 11/15/2026.
rungs = [0.10 0.50 0.90 0.99 0.999 0.9999];
pf_pinned(1,{'DD'},{'none','box','box0','faces','faces0'},[0 1 2],rungs,'c0_dd');
pf_pinned(2,{'DD','DD'},{'none','box','box0','product','product0','faces','faces0','facesS','corners','corners0'},[0 1],rungs,'c1_dddd');
disp('## C01 DONE');
