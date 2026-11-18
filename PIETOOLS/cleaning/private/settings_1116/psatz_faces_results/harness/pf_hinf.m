function R = pf_hinf(PIEfun,gex,cfgs,dws,storages,tier,tag)
% R = PF_HINF(PIEFUN,GEX,CFGS,DWS,STORAGES,TIER,TAG) the coercive H-infinity
% gain executive (PIETOOLS_Hinf_gain_sop, form 'P') on the PIE returned by
% PIEFUN() (a function handle, built after the path is set), closed-form
% gain GEX, for every term-set configuration in CFGS (PF_CFG), every fixed
% slack weight raise dw in DWS and every storage in STORAGES ('stock',
% 'def'). The term set enters through the tier knob lpi.ineq. Records
% gamma, the gap gamma/gex - 1 and the executive's certified flag; results
% appended to results/hinf_<TAG>.mat. SCRATCH, MMP, 11/15/2026.
R = struct('storage',{},'cfg',{},'o',{},'pd',{},'dw',{},'status',{},'certified',{},'gam',{},'gap',{}, ...
           'ndv',{},'m',{},'t',{},'err',{});
for is = 1:numel(storages)
    S = pf_setup(strcmp(storages{is},'def'));
    out = sprintf('%s/results/hinf_%s.mat',S,tag);
    PIE = PIEfun();
    for ic = 1:numel(cfgs)
        c = cfgs(ic);   o = NaN;    if ~isempty(c.offs),  o = c.offs(1);    end
        for dw = dws
            r = struct('storage',storages{is},'cfg',c.name,'o',o,'pd',c.pd,'dw',dw,'status','','certified',false, ...
                       'gam',NaN,'gap',NaN,'ndv',NaN,'m',NaN,'t',NaN,'err','');
            t0 = tic;
            try
                st = lpisettings_sop(tier);
                st.sos_opts.solver = 'mosek';   st.sop.verbose = false;
                st.sop.dw_P = dw;   st.sop.dw_Q = dw;   st.sop.max_raise = [0 0 0];
                st.lpi.ineq.psatz = 0;  st.lpi.ineq.psatz_offset = 0;
                st.lpi.ineq.psatz_rows = c.rows;    st.lpi.ineq.psatz_rows_offset = c.offs;
                st.lpi.ineq.psatz_rows_perdir = c.pd;
                [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,st,'P');
                r.t = toc(t0);
                r.status = info.status;     r.certified = logical(info.certified);
                r.gam = gam;    r.gap = gam/gex - 1;
                r.ndv = info.shape.ndv;     r.m = info.shape.m;
            catch ME
                r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
            end
            R(end+1) = r;                                                   %#ok<AGROW>
            fprintf('@@ HINF storage=%s cfg=%s o=%g pd=%d dw=%d %s cert=%d gam=%.7f gap=%+.2e ndv=%d m=%d t=%.0f %s\n', ...
                    storages{is},c.name,o,c.pd,dw,r.status,r.certified,r.gam,r.gap,r.ndv,r.m,r.t,r.err);
            save(out,'R');
        end
    end
end
fprintf('## HINF DONE %s\n',tag);
end
