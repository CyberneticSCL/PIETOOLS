% PF_RUN_C2F  The half step where the theory is exact, and on the pinned
% plant. (a) 1-D heat DD (kappa* = pi^2), Q at (2,wQ): the pair (box o1)
% at wQ = 1 and 2, box0 and faces0 at wQ = 1, rtol 1e-5: Markov-Lukacs
% says the pair is complete at even degree 2w and the faces at odd degree
% 2w + 1, so faces0 at w should lie between the pair at w and at w + 1.
% (b) the pinned-P 2-D DD x DD plant at dw 0 (w [1 1]) on finer rungs
% between 0.50 and 0.90 for every set, to place the half step there too.
% SCRATCH, MMP, 11/15/2026.
bo = struct('rtol',1e-5,'maxsolve',30,'tmax',600);
C = pf_cfg({'box','box0','faces0'},1);
pf_heat(1,{'DD'},0,0.1,C,[1 1],1,'c2f_1d_wq1',bo);
C2 = pf_cfg({'box'},1);
pf_heat(1,{'DD'},0,0.1,C2,[1 1],2,'c2f_1d_wq2',bo);
rungs = [0.55 0.60 0.65 0.70 0.75 0.80 0.85];
pf_pinned(2,{'DD','DD'},{'box','box0','product0','faces','faces0','facesS','corners','corners0'},0,rungs,'c2f_pinned_fine');
disp('## C2F DONE');
