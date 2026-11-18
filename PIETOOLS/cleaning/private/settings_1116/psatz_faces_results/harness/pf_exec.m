function R = pf_exec(sysname,pf,pvals,cfgs,dws,storages,tier,tag,form)
% R = PF_EXEC(SYS,PF,PVALS,CFGS,DWS,STORAGES,TIER,TAG) the executive path
% (PIETOOLS_stability_sop, form 'P') on the mixed-domain system SYS
% (MIXED_DOMAIN_SYSTEMS) at parameter PF = each of PVALS, for every term-set
% configuration in CFGS (PF_CFG), every fixed slack weight raise dw in DWS
% (no raise loop: dw_feas = dw, max_raise_feas = 0) and every storage in
% STORAGES ('stock': the repo's settings_2d storage; 'def': the route
% shadow, POSCOPVAR with DEG omitted). The term set enters through the tier
% knob lpi.ineq (psatz 0 + rows), read by the scratch lpi_ineq_sop. The
% verdict is the executive's own (info.certified). Results appended to
% results/exec_<TAG>.mat after every run. SCRATCH, MMP, 11/15/2026.
if nargin<9 || isempty(form),   form = 'P';     end
R = struct('sys',{},'pf',{},'p',{},'storage',{},'cfg',{},'o',{},'pd',{},'dw',{},'status',{},'certified',{}, ...
           'ndv',{},'m',{},'t',{},'err',{});
for is = 1:numel(storages)
    S = pf_setup(strcmp(storages{is},'def'));
    out = sprintf('%s/results/exec_%s.mat',S,tag);
    Sys = mixed_domain_systems(sysname,struct(pf,pvals(1)));                %#ok<NASGU>
    for p = pvals
        Sys = mixed_domain_systems(sysname,struct(pf,p));
        PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
        for ic = 1:numel(cfgs)
            c = cfgs(ic);   o = NaN;    if ~isempty(c.offs),  o = c.offs(1);    end
            for dw = dws
                r = struct('sys',sysname,'pf',pf,'p',p,'storage',storages{is},'cfg',c.name,'o',o,'pd',c.pd,'dw',dw, ...
                           'status','','certified',false,'ndv',NaN,'m',NaN,'t',NaN,'err','');
                try
                    st = lpisettings_sop(tier);
                    st.sos_opts.solver = 'mosek';   st.sop.verbose = false;
                    st.sop.dw_feas = dw*[1 1 1];    st.sop.max_raise_feas = [0 0 0];
                    st.lpi.ineq.psatz = 0;  st.lpi.ineq.psatz_offset = 0;
                    st.lpi.ineq.psatz_rows = c.rows;    st.lpi.ineq.psatz_rows_offset = c.offs;
                    st.lpi.ineq.psatz_rows_perdir = c.pd;
                    t0 = tic;
                    [~,~,info] = PIETOOLS_stability_sop(PIE,st,form);
                    r.t = toc(t0);
                    r.status = info.status;     r.certified = logical(info.certified);
                    r.ndv = info.shape.ndv;     r.m = info.shape.m;
                catch ME
                    r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
                end
                R(end+1) = r;                                               %#ok<AGROW>
                fprintf('@@ EXEC %s %s=%.6f storage=%s cfg=%s o=%g pd=%d dw=%d %s cert=%d ndv=%d m=%d t=%.0f %s\n', ...
                        sysname,pf,p,storages{is},c.name,o,c.pd,dw,r.status,r.certified,r.ndv,r.m,r.t,r.err);
                save(out,'R');
            end
        end
    end
end
fprintf('## EXEC DONE %s\n',tag);
end
