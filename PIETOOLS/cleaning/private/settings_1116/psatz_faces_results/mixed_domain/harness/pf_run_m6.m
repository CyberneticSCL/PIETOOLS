% PF_RUN_M6  S18 (tau = 0.5 tau*) with the delay direction s2 split (lower
% and upper kernels, 'sep10') and the storage degrees raised above every
% tier's 2-D entry (mult 1, int 2 per variable, joint 4, for every space;
% PF_POS override read by PF_POS_SEP), facesS, dw 0 and 1, MOSEK; then both
% directions split ('sep00') at dw 1. SCRATCH, MMP, 11/16/2026.
setappdata(0,'PF_POS',struct('mult',[1 1 1],'int',[2 2 2],'joint',[4 4 4],'sep',0));
C = pf_cfg({'facesS'},2);
tau43 = 1.034726;
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,[0 1],{'sep10'},'light','m6_s18_deg2');
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,1,{'sep00'},'light','m6_s18_deg2_sep00');
rmappdata(0,'PF_POS');
disp('## M6 DONE');
