function Out = pf_lib(ids,cfgs,storages,tier,tag)
% OUT = PF_LIB(IDS,CFGS,STORAGES,TIER,TAG) the 2-D rows of the PDE example
% library (catalog of the route4 library harness) through the container
% executives with the raise loop as the tier sets it: stability 'P' on
% 33, 34, 35, 37, 39, 40; H-infinity gain 'P' on 45, 46, 47; the H-infinity
% estimator on 48; for every term-set configuration in CFGS (tier knob
% lpi.ineq, scratch lpi_ineq_sop) and storage ('stock' / 'def'). Results
% appended to results/lib_<TAG>.mat. SCRATCH, MMP, 11/15/2026.
H = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad';
LB = [H '/psatz_mix/route4/lb'];
Out = struct('id',{},'task',{},'storage',{},'cfg',{},'o',{},'pd',{},'status',{},'certified',{},'dw',{},'gam',{},'ndv',{},'m',{},'t',{},'err',{});
for is = 1:numel(storages)
    S = pf_setup(strcmp(storages{is},'def'));
    addpath(LB);
    out = sprintf('%s/results/lib_%s.mat',S,tag);
    L = load([LB '/catalog.mat']);  C = L.C;
    rows = C(strcmp({C.lib},'PDE') & ismember({C.id},ids));
    [~,iu] = unique({rows.id},'stable');    rows = rows(iu);
    for i = 1:numel(rows)
        mrow = rows(i);
        try
            params = {};
            if mrow.idxarg,   PDE_t = feval(mrow.fn,str2double(mrow.id),0,params);
            else,             PDE_t = feval(mrow.fn,0,params);
            end
            PDE = initialize_PIETOOLS_PDE(PDE_t,true);
            PIE = initialize(convert(PDE,{'pie','silent'}));
        catch ME
            fprintf('@@ LIB %s build ERROR %s\n',mrow.id,strrep(ME.message,newline,' '));   continue
        end
        switch mrow.id
            case {'33','34','35','37','39','40'},   task = 'stabP';
            case {'45','46','47'},                  task = 'hinfP';
            case '48',                              task = 'est';
            otherwise,                              task = 'stabP';
        end
        for ic = 1:numel(cfgs)
            c = cfgs(ic);   o = NaN;    if ~isempty(c.offs),  o = c.offs(1);    end
            r = struct('id',mrow.id,'task',task,'storage',storages{is},'cfg',c.name,'o',o,'pd',c.pd, ...
                       'status','','certified',false,'dw',NaN,'gam',NaN,'ndv',NaN,'m',NaN,'t',NaN,'err','');
            t0 = tic;
            try
                st = lpisettings_sop(tier);
                st.sos_opts.solver = 'mosek';   st.sop.verbose = false;
                st.lpi.ineq.psatz = 0;  st.lpi.ineq.psatz_offset = 0;
                st.lpi.ineq.psatz_rows = c.rows;    st.lpi.ineq.psatz_rows_offset = c.offs;
                st.lpi.ineq.psatz_rows_perdir = c.pd;
                gam = NaN;
                switch task
                    case 'stabP',   [~,~,info] = PIETOOLS_stability_sop(PIE,st,'P');
                    case 'hinfP',   [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,st,'P');
                    case 'est',     [~,~,gam,~,~,info] = PIETOOLS_Hinf_estimator_sop(PIE,st);
                end
                r.t = toc(t0);  r.status = info.status;     r.certified = logical(info.certified);
                r.dw = info.dw; r.gam = gam;    r.ndv = info.shape.ndv;     r.m = info.shape.m;
            catch ME
                r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
            end
            Out(end+1) = r;                                                 %#ok<AGROW>
            fprintf('@@ LIB %-3s %-6s storage=%s cfg=%-9s o=%g pd=%d %-10s cert=%d dw=%g gam=%.10g ndv=%d m=%d t=%.0f %s\n', ...
                    mrow.id,task,storages{is},c.name,o,c.pd,r.status,r.certified,r.dw,r.gam,r.ndv,r.m,r.t,r.err);
            save(out,'Out');
        end
    end
end
fprintf('## LIB DONE %s\n',tag);
end
