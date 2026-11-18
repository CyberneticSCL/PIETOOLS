% PF_RUN_M9  Attribution of the S18 certificate to the storage degree: the
% tier's own DEG-omitted storage ('def': the pure multiplier and the pure
% full integral of the 2-D space, separable) with lpi.pos raised to mult 1,
% int 2, joint 4 (PF_POS, read by PF_EXEC_M into st.lpi.pos), facesS, dw 0,
% at 0.5 tau*; then the same at 1.1 tau* (must fail). SCRATCH, MMP, 11/16/2026.
setappdata(0,'PF_POS',struct('mult',[1 1 1],'int',[2 2 2],'joint',[4 4 4],'sep',0));
C = pf_cfg({'facesS'},2);
tau43 = 1.034726;
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1),struct('tau',1.1*tau43,'sgn',-1)},C,0,{'def'},'light','m9_s18_def_deg2');
rmappdata(0,'PF_POS');
disp('## M9 DONE');
