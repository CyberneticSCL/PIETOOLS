function R = pf_dmax(sysname,plist,cfgs,dws,storages,tier,tag,dcap)
% R = PF_DMAX(SYS,PLIST,CFGS,DWS,STORAGES,TIER,TAG,DCAP) the decay-margin
% program of the stability LPI on the executive path (the storage, slack
% sizing and caps of PIETOOLS_stability_sop, form 'P'): with P = Pm + I (Pm
% the storage, so P >= I), maximize delta subject to
%       -(T'*P*A + A'*P*T) - delta*(T'*T) >= 0   (LPI_INEQ_SOP like P, dw pinned)
%       0 <= delta <= DCAP (default 1).
% By homogeneity in P, a strict certificate at these degrees gives delta* =
% DCAP; 'optimal' with delta* = 0 means that no coercive P of this storage
% has a derivative operator with a margin (a program without interior, where
% the executive's feasibility form returns UNKNOWN); 'infeasible' that no
% P >= I has a nonpositive derivative at all. Storages, term sets and the
% solver flags as in PF_EXEC_M. Results appended to results/dmax_<TAG>.mat.
% SCRATCH, MMP, 11/16/2026.
if nargin<8 || isempty(dcap),   dcap = 1;   end
R = struct('sys',{},'par',{},'storage',{},'cfg',{},'dw',{},'status',{},'delta',{},'ndv',{},'t',{}, ...
           'numerr',{},'feasratio',{},'relres',{},'eta',{},'err',{});
for is = 1:numel(storages)
    sn = storages{is};
    if strcmp(sn,'stock'),          S = pf_setup(false);
    elseif strcmp(sn,'def'),        S = pf_setup('route');
    elseif startsWith(sn,'sep'),    S = pf_setup('route_sep');     setappdata(0,'PF_SEP',sn(4:end)=='1');
    else,   error('pf_dmax:storage','unknown storage %s',sn);
    end
    out = sprintf('%s/results/dmax_%s.mat',S,tag);
    for ip = 1:numel(plist)
        p = plist{ip};
        Sys = mixed_domain_systems(sysname,p);
        PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
        PIE = exec_pie_sop(PIE);
        assert(~(isa(PIE.K,'copvar') && size(PIE.K,1)>0),'pf_dmax: a constraint K x = 0 is not handled here');
        Tm = PIE.T;     Am = PIE.A;
        for ic = 1:numel(cfgs)
            c = cfgs(ic);
            for dw = dws
                r = struct('sys',sysname,'par',p,'storage',sn,'cfg',c.name,'dw',dw,'status','','delta',NaN, ...
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
                    [prog,dv] = lpidecvar(prog,'delta');
                    prog = lpi_ineq_sop(prog,dv);   prog = lpi_ineq_sop(prog,dcap-dv);
                    prog = lpisetobj(prog,-dv);
                    [prog,Pm] = tl.storage(prog,Tm,'out');     Pm = Pm + tl.eye(1,1);
                    Dm = Tm'*Pm*Am + Am'*Pm*Tm + dv*(Tm'*Tm);
                    o = struct('like',Pm,'dw',dw,'margin',st.sop.margin,'cap',tl.capfor(Dm));
                    [prog,~,~] = lpi_ineq_sop(prog,-Dm,o);
                    r.ndv = size(prog.decvartable,1);
                    prog = lpisolve(prog,st.sos_opts);
                    cl = lpi_classify_sop(prog,dv);
                    r.t = toc(t0);
                    r.status = cl.status;   r.delta = cl.obj;   r.numerr = cl.numerr;
                    r.feasratio = cl.feasratio;     r.relres = cl.relres;  r.eta = cl.eta;
                catch ME
                    r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
                end
                R(end+1) = r;                                               %#ok<AGROW>
                fprintf('@@ DMAX %s par=%s storage=%s cfg=%s dw=%d %s delta=%.6g ndv=%d t=%.0f numerr=%g feas=%.3g relres=%.1e eta=%.1e %s\n', ...
                        sysname,jsonencode(p),sn,c.name,dw,r.status,r.delta,r.ndv,r.t,r.numerr,r.feasratio,r.relres,r.eta,r.err);
                save(out,'R');
            end
        end
    end
end
fprintf('## DMAX DONE %s\n',tag);
end
