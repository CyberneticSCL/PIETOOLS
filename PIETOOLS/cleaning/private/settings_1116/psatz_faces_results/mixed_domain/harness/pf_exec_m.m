function R = pf_exec_m(sysname,plist,cfgs,dws,storages,tier,tag,solver,form)
% R = PF_EXEC_M(SYS,PLIST,CFGS,DWS,STORAGES,TIER,TAG,SOLVER,FORM) is PF_EXEC
% for the S17/S18 diagnostics: PLIST a cell of parameter structs of the
% system SYS (MIXED_DOMAIN_SYSTEMS); STORAGES a cell of names, 'stock' (the
% repo's settings_2d storage), 'def' (the route shadow: POSCOPVAR with DEG
% omitted, the N-D default) or 'sepAB' (the route_sep shadow: every basis
% operator at the light degrees, registry direction d split into lower and
% upper integrals where digit d is 0 and separable (full integral) where it
% is 1, PF_POS_SEP); SOLVER 'mosek' (default) or 'sedumi'. The term set
% enters through the tier knob lpi.ineq (psatz 0 + rows); dw is pinned (no
% raise loop). The verdict is the executive's (info.certified); the solver
% flags numerr, feasratio, relres, eta are read from prog.solinfo.info.
% Results appended to results/exec_<TAG>.mat after every run.
% SCRATCH, MMP, 11/16/2026.
if nargin<8 || isempty(solver),     solver = 'mosek';   end
if nargin<9 || isempty(form),       form = 'P';         end
R = struct('sys',{},'par',{},'storage',{},'cfg',{},'o',{},'pd',{},'dw',{},'status',{},'certified',{}, ...
           'ndv',{},'m',{},'t',{},'solver',{},'numerr',{},'feasratio',{},'relres',{},'eta',{},'err',{});
getf = @(s,f) iff(isstruct(s) && isfield(s,f) && ~isempty(s.(f)),s,f);
for is = 1:numel(storages)
    sn = storages{is};
    if strcmp(sn,'stock'),          S = pf_setup(false,strcmp(form,'Q'));
    elseif strcmp(sn,'def'),        S = pf_setup('route',strcmp(form,'Q'));
    elseif startsWith(sn,'sep'),    S = pf_setup('route_sep',strcmp(form,'Q'));     setappdata(0,'PF_SEP',sn(4:end)=='1');
    elseif strcmp(sn,'row'),        S = pf_setup('route_row',strcmp(form,'Q'));     % the row in appdata PF_ROW
    else,   error('pf_exec_m:storage','unknown storage %s',sn);
    end
    out = sprintf('%s/results/exec_%s.mat',S,tag);
    for ip = 1:numel(plist)
        p = plist{ip};
        Sys = mixed_domain_systems(sysname,p);
        PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
        for ic = 1:numel(cfgs)
            c = cfgs(ic);   o = NaN;    if ~isempty(c.offs),  o = c.offs(1);    end
            for dw = dws
                r = struct('sys',sysname,'par',p,'storage',sn,'cfg',c.name,'o',o,'pd',c.pd,'dw',dw, ...
                           'status','','certified',false,'ndv',NaN,'m',NaN,'t',NaN,'solver',solver, ...
                           'numerr',NaN,'feasratio',NaN,'relres',NaN,'eta',NaN,'err','');
                t0 = tic;
                try
                    st = lpisettings_sop(tier);
                    if isappdata(0,'PF_POS') && ~isempty(getappdata(0,'PF_POS')),  st.lpi.pos = getappdata(0,'PF_POS');  end    % storage degree override (M9)
                    st.sos_opts.solver = solver;    st.sop.verbose = false;
                    if isappdata(0,'PF_EPNEG') && ~isempty(getappdata(0,'PF_EPNEG')),  st.epneg = getappdata(0,'PF_EPNEG');  end    % decay rate override (M10, the FP form)
                    st.sop.dw_feas = dw*[1 1 1];    st.sop.max_raise_feas = [0 0 0];
                    st.lpi.ineq.psatz = 0;  st.lpi.ineq.psatz_offset = 0;
                    st.lpi.ineq.psatz_rows = c.rows;    st.lpi.ineq.psatz_rows_offset = c.offs;
                    st.lpi.ineq.psatz_rows_perdir = c.pd;
                    [prog,~,info] = PIETOOLS_stability_sop(PIE,st,form);
                    r.t = toc(t0);
                    r.status = info.status;     r.certified = logical(info.certified);
                    r.ndv = info.shape.ndv;     r.m = info.shape.m;
                    si = [];
                    if isfield(prog,'solinfo') && isfield(prog.solinfo,'info'),    si = prog.solinfo.info;     end
                    r.numerr = getf(si,'numerr');   r.feasratio = getf(si,'feasratio');
                    r.relres = getf(si,'relres');   r.eta = getf(si,'eta');
                catch ME
                    r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
                end
                R(end+1) = r;                                               %#ok<AGROW>
                fprintf('@@ EXECM %s par=%s storage=%s cfg=%s o=%g pd=%d dw=%d %s cert=%d ndv=%d m=%d t=%.0f numerr=%g feas=%.3g relres=%.1e eta=%.1e %s %s\n', ...
                        sysname,jsonencode(p),sn,c.name,o,c.pd,dw,r.status,r.certified,r.ndv,r.m,r.t, ...
                        r.numerr,r.feasratio,r.relres,r.eta,solver,r.err);
                save(out,'R');
            end
        end
    end
end
fprintf('## EXECM DONE %s\n',tag);
end

function v = iff(c,s,f)
% V = IFF(C,S,F) the field F of S when C, else NaN.
if c,   v = s.(f);  else,   v = NaN;    end
end
