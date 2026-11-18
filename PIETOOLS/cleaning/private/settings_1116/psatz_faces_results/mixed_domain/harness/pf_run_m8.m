% PF_RUN_M8  Attribution and soundness of the S18 certificate of M6 (s2
% split, raised degrees mult 1 / int 2 / joint 4, facesS, dw 0): the
% separable storage at the same degrees ('sep11', control for the split) at
% 0.5 tau*; the certifying storage at 1.1 tau* (unstable: must not certify)
% and at 0.9 tau*. SCRATCH, MMP, 11/16/2026.
setappdata(0,'PF_POS',struct('mult',[1 1 1],'int',[2 2 2],'joint',[4 4 4],'sep',0));
C = pf_cfg({'facesS'},2);
tau43 = 1.034726;
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,0,{'sep11'},'light','m8_s18_sep11_deg2');
pf_exec_m('S18_lib43_interior_delay',{struct('tau',1.1*tau43,'sgn',-1),struct('tau',0.9*tau43,'sgn',-1)},C,0,{'sep10'},'light','m8_s18_sep10_deg2_margin');
rmappdata(0,'PF_POS');
disp('## M8 DONE');
