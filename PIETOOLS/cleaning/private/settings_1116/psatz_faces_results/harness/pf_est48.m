function Out = pf_est48(cfgs,dws,storages,tier,tag)
% OUT = PF_EST48(CFGS,DWS,STORAGES,TIER,TAG) the H-infinity estimator of
% library PDE 48 (2-D) with the slack weight pinned (dw_P = dw_Q = dw,
% max_raise 0), every term-set configuration in CFGS (tier knob lpi.ineq),
% both storages: the estimator gain bound, the executive's verdict, sizes
% and times. The only 2-D library row on which the term sets differed at
% light (C6: box unknown, faces 0.00465, facesS 0.00618, product 0.00646).
% Results appended to results/est48_<TAG>.mat. SCRATCH, MMP, 11/15/2026.
H = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad';
LB = [H '/psatz_mix/route4/lb'];
Out = struct('storage',{},'cfg',{},'o',{},'pd',{},'dw',{},'status',{},'certified',{},'gam',{},'ndv',{},'m',{},'t',{},'err',{});
for is = 1:numel(storages)
    S = pf_setup(strcmp(storages{is},'def'));
    addpath(LB);
    out = sprintf('%s/results/est48_%s.mat',S,tag);
    L = load([LB '/catalog.mat']);  C = L.C;
    mrow = C(strcmp({C.lib},'PDE') & strcmp({C.id},'48'));    mrow = mrow(1);
    params = {};
    if mrow.idxarg,   PDE_t = feval(mrow.fn,48,0,params);  else,   PDE_t = feval(mrow.fn,0,params);   end
    PDE = initialize_PIETOOLS_PDE(PDE_t,true);
    PIE = initialize(convert(PDE,{'pie','silent'}));
    for ic = 1:numel(cfgs)
        c = cfgs(ic);   o = NaN;    if ~isempty(c.offs),  o = c.offs(1);    end
        for dw = dws
            r = struct('storage',storages{is},'cfg',c.name,'o',o,'pd',c.pd,'dw',dw,'status','','certified',false, ...
                       'gam',NaN,'ndv',NaN,'m',NaN,'t',NaN,'err','');
            t0 = tic;
            try
                st = lpisettings_sop(tier);
                st.sos_opts.solver = 'mosek';   st.sop.verbose = false;
                st.sop.dw_P = dw;   st.sop.dw_Q = dw;   st.sop.max_raise = [0 0 0];
                st.lpi.ineq.psatz = 0;  st.lpi.ineq.psatz_offset = 0;
                st.lpi.ineq.psatz_rows = c.rows;    st.lpi.ineq.psatz_rows_offset = c.offs;
                st.lpi.ineq.psatz_rows_perdir = c.pd;
                [~,~,gam,~,~,info] = PIETOOLS_Hinf_estimator_sop(PIE,st);
                r.t = toc(t0);  r.status = info.status;     r.certified = logical(info.certified);
                r.gam = gam;    r.ndv = info.shape.ndv;     r.m = info.shape.m;
            catch ME
                r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
            end
            Out(end+1) = r;                                                 %#ok<AGROW>
            fprintf('@@ EST48 storage=%s cfg=%-9s o=%g pd=%d dw=%d %-10s cert=%d gam=%.10g ndv=%d m=%d t=%.0f %s\n', ...
                    storages{is},c.name,o,c.pd,dw,r.status,r.certified,r.gam,r.ndv,r.m,r.t,r.err);
            save(out,'Out');
        end
    end
end
fprintf('## EST48 DONE %s\n',tag);
end
