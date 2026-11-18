% PF_RUN_C2B  Phase 5 (ladder): the Putinar box set one weight degree up
% on the 2-D heat benchmark, Q at (2,3), to test whether it reaches what
% the faces reach at (2,2) (T6, H1). Budget 12 solves / 2400 s per
% configuration. SCRATCH, MMP, 11/15/2026.
C = pf_cfg({'box'},2);
pf_heat(2,{'DD','DN'},0,0.1,C,[1 1],3,'c2b_wq3',struct('maxsolve',12,'tmax',2400));
disp('## C2B DONE');
