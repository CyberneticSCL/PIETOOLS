% PF_RUN_FIXED  Positivity certification of fixed operators on the mixed
% spaces of S17, S3 and S18 (PF_FIXED): the PDE-norm form T'*T, its
% restrictions to the 1-D space (T'*E1*T) and to the 2-D space (T'*E2*T),
% and A'*A; for S18 also the certified storage P* of M8 (every separable
% basis operator at the raised degrees) and the certified derivative
% operator -(T'*P*A + A'*P*T), as fixed targets. Options: the reader's
% defaults; no terms; raised dD, dw. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
tau43 = 1.034726;
sys = {'S17_transport_edge_diffusion',struct('r2',0.5*pi^2); 'S3_edge1d_heat2d',struct('r',0.9*pi^2); 'S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1)};
ops = struct('name',{},'X',{});
for i = 1:size(sys,1)
    Sys = mixed_domain_systems(sys{i,1},sys{i,2});
    PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
    T = PIE.T;  A = PIE.A;
    st = exec_sop_settings(lpisettings_sop('light'));   tl = exec_tools_sop(PIE,st);
    E1 = tl.eye(1,[0 1 0 0]);   E2 = tl.eye(1,[0 0 0 1]);   % weights by variable subset: L2[s1] entry 2, L2[s1,s2] entry 4
    nm = sys{i,1}(1:3);
    ops(end+1) = struct('name',[nm '_TtT'],'X',T'*T);                       %#ok<AGROW>
    ops(end+1) = struct('name',[nm '_TtE1T'],'X',T'*E1*T);                  %#ok<AGROW>
    ops(end+1) = struct('name',[nm '_TtE2T'],'X',T'*E2*T);                  %#ok<AGROW>
    ops(end+1) = struct('name',[nm '_AtA'],'X',A'*A);                       %#ok<AGROW>
    if strncmp(nm,'S18',3)
        % the certified storage of M8 (sep11, raised degrees) and its derivative operator
        pf_setup('route_sep');  setappdata(0,'PF_SEP',[true true]);
        setappdata(0,'PF_POS',struct('mult',[1 1 1],'int',[2 2 2],'joint',[4 4 4],'sep',0));
        st2 = lpisettings_sop('light');     st2.sos_opts.solver = 'mosek';  st2.sop.verbose = false;
        st2.sop.dw_feas = [0 0 0];  st2.sop.max_raise_feas = [0 0 0];
        st2.lpi.pos = getappdata(0,'PF_POS');
        C = pf_cfg({'facesS'},2);   c = C(1);
        st2.lpi.ineq.psatz = 0; st2.lpi.ineq.psatz_offset = 0;  st2.lpi.ineq.psatz_rows = c.rows;
        st2.lpi.ineq.psatz_rows_offset = c.offs;    st2.lpi.ineq.psatz_rows_perdir = c.pd;
        [~,Pst,info] = PIETOOLS_stability_sop(PIE,st2,'P');
        fprintf('@@ FIX S18 storage solve: %s cert=%d ndv=%d\n',info.status,info.certified,info.shape.ndv);
        save(sprintf('%s/results/s18_Pstar.mat',S),'Pst','info');
        rmappdata(0,'PF_POS');  pf_setup(false);
        ops(end+1) = struct('name','S18_Pstar','X',Pst);                    %#ok<AGROW>
        ops(end+1) = struct('name','S18_negD','X',-(T'*Pst*A + A'*Pst*T));  %#ok<AGROW>
    end
end
optl = {struct(), struct('psatz','none'), struct('dD',2,'dw',2)};
pf_fixed('mixed',ops,optl,'light');
disp('## FIXED DONE');
