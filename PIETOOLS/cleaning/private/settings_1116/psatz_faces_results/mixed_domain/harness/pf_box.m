function R = pf_box(sysname,pdef,cfgs,dw,storages,tier,tag,epl)
% R = PF_BOX(SYS,PDEF,CFGS,DW,STORAGES,TIER,TAG,EPL) the stability LPI of the
% executive path with the positive operator bounded on both sides,
%       ep I <= P <= I,    -(T'*P*A + A'*P*T) >= 0,
% for each ep in EPL: P = Ps + ep I with Ps the storage, the upper bound as
% the inequality (1-ep) I - Ps >= 0 (LPI_INEQ_SOP sized like Ps), the
% derivative inequality sized like P at the pinned dw. The bounds fix the
% scale of P, so a program that is infeasible only through states with a
% boundary layer of width ep/||P11|| (S17, D4: S17_obstruction.md) becomes
% strongly infeasible at ep of order 1 and the solver returns a Farkas
% certificate; a genuine certificate with condition number <= 1/ep stays
% feasible. SYS: a name of MIXED_DOMAIN_SYSTEMS with PDEF its parameter
% struct, or 'D4' (PDEF the reaction r2). Results appended to
% results/box_<TAG>.mat. SCRATCH, MMP, 11/16/2026.
R = struct('sys',{},'storage',{},'cfg',{},'dw',{},'ep',{},'status',{},'ndv',{},'t',{}, ...
           'numerr',{},'feasratio',{},'relres',{},'eta',{},'err',{});
for is = 1:numel(storages)
    sn = storages{is};
    if strcmp(sn,'stock'),          S = pf_setup(false);
    elseif strcmp(sn,'def'),        S = pf_setup('route');
    elseif startsWith(sn,'sep'),    S = pf_setup('route_sep');     setappdata(0,'PF_SEP',sn(4:end)=='1');
    else,   error('pf_box:storage','unknown storage %s',sn);
    end
    out = sprintf('%s/results/box_%s.mat',S,tag);
    if strcmp(sysname,'D4')
        pvar s1 s2 t
        clear stateNameGenerator
        x1 = pde_var(1,s1,[0,1]);   x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
        PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s2,2)+pdef*x2; ...
               subs(x1,s1,0)==0; subs(x1,s1,1)==0; subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
    else
        Sys = mixed_domain_systems(sysname,pdef);   PDE = Sys.PDE;
    end
    PIE = initialize(convert(PDE,{'pie','silent'}));
    PIE = exec_pie_sop(PIE);
    assert(~(isa(PIE.K,'copvar') && size(PIE.K,1)>0),'pf_box: a constraint K x = 0 is not handled here');
    Tm = PIE.T;     Am = PIE.A;
    for ic = 1:numel(cfgs)
        c = cfgs(ic);
        for ep = epl
            r = struct('sys',sysname,'storage',sn,'cfg',c.name,'dw',dw,'ep',ep,'status','', ...
                       'ndv',NaN,'t',NaN,'numerr',NaN,'feasratio',NaN,'relres',NaN,'eta',NaN,'err','');
            t0 = tic;
            try
                st = lpisettings_sop(tier);
                st.sos_opts.solver = 'mosek';   st.sop.verbose = false;
                st.sop.dw_feas = dw*[1 1 1];    st.sop.max_raise_feas = [0 0 0];
                st.lpi.ineq.psatz = 0;  st.lpi.ineq.psatz_offset = 0;
                st.lpi.ineq.psatz_rows = c.rows;    st.lpi.ineq.psatz_rows_offset = c.offs;
                st.lpi.ineq.psatz_rows_perdir = c.pd;
                st = exec_sop_settings(st);
                tl = exec_tools_sop(PIE,st);    st.sop = tl.sop;
                prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
                [prog,Ps] = tl.storage(prog,Tm,'out');
                Pm = Ps + tl.eye(1,ep);
                ou = struct('like',Ps,'dw',dw,'margin',st.sop.margin);
                [prog,~,~] = lpi_ineq_sop(prog,tl.eye(1,1-ep) - Ps,ou);
                Dm = Tm'*Pm*Am + Am'*Pm*Tm;
                o = struct('like',Pm,'dw',dw,'margin',st.sop.margin,'cap',tl.capfor(Dm));
                [prog,~,~] = lpi_ineq_sop(prog,-Dm,o);
                r.ndv = size(prog.decvartable,1);
                prog = lpisolve(prog,st.sos_opts);
                cl = lpi_classify_sop(prog,[]);
                r.t = toc(t0);
                r.status = cl.status;   r.numerr = cl.numerr;
                r.feasratio = cl.feasratio;     r.relres = cl.relres;  r.eta = cl.eta;
            catch ME
                r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
            end
            R(end+1) = r;                                                   %#ok<AGROW>
            fprintf('@@ BOX %s storage=%s cfg=%s dw=%d ep=%g %s ndv=%d t=%.0f numerr=%g feas=%.3g relres=%.1e eta=%.1e %s\n', ...
                    sysname,sn,c.name,dw,ep,r.status,r.ndv,r.t,r.numerr,r.feasratio,r.relres,r.eta,r.err);
            save(out,'R');
        end
    end
end
fprintf('## BOX DONE %s\n',tag);
end
