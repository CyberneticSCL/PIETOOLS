function r = pf_exec_row(label,PIE,task,storage,tier,extra)
% R = PF_EXEC_ROW(LABEL,PIE,TASK,STORAGE,TIER,EXTRA) one executive call on
% the container PIE under a storage: TASK 'stabP' (PIETOOLS_stability_sop
% 'P'), 'hinfP' (PIETOOLS_Hinf_gain_sop 'P'), 'hinfQ', 'est'
% (PIETOOLS_Hinf_estimator_sop); STORAGE 'stock' (the repo's), 'def' (the
% N-D default route shadow) or a row name of PF_ROWS (the row shadow);
% TIER the lpisettings_sop tier; EXTRA (optional) a struct of settings
% overrides (fields applied to st.sop, st.lpi.ineq or st). The executive's
% own raise loop and terms. Prints '@@ EX' and returns the record.
% SCRATCH, MMP, 11/16/2026.
if nargin<6,    extra = struct();   end
switch storage
    case 'stock',   pf_setup(false);
    case 'def',     pf_setup('route');
    otherwise,      pf_setup('route_row');  setappdata(0,'PF_ROW',pf_rows(storage));
end
r = struct('label',label,'task',task,'storage',storage,'tier',tier,'status','','certified',false,'dw',NaN,'gam',NaN, ...
           'ndv',NaN,'m',NaN,'t',NaN,'msg','');
t0 = tic;
try
    st = lpisettings_sop(tier);     st.sos_opts.solver = 'mosek';   st.sop.verbose = false;
    fn = fieldnames(extra);
    for i = 1:numel(fn)
        f = fn{i};
        if isfield(st.sop,f),           st.sop.(f) = extra.(f);
        elseif isfield(st.lpi.ineq,f),  st.lpi.ineq.(f) = extra.(f);
        else,                           st.(f) = extra.(f);
        end
    end
    gam = NaN;
    switch task
        case 'stabP',   [~,~,info] = PIETOOLS_stability_sop(PIE,st,'P');
        case 'hinfP',   [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,st,'P');
        case 'hinfQ',   [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,st,'Q');
        case 'est',     [txt,~,~,gam,~,~,info] = evalc('PIETOOLS_Hinf_estimator_sop(PIE,st)');  %#ok<ASGLU>
        otherwise,      error('pf_exec_row:task','unknown task %s',task);
    end
    r.status = info.status;     r.certified = logical(info.certified);   r.dw = info.dw;
    r.gam = gam;    r.ndv = info.shape.ndv;     r.m = info.shape.m;
catch ME
    r.status = 'ERROR';     r.msg = strrep(ME.message,newline,' ');
end
r.t = toc(t0);
fprintf('@@ EX %s task=%s storage=%s tier=%s %s cert=%d dw=%g gam=%.9g ndv=%d m=%d t=%.0f %s\n', ...
        label,task,storage,tier,r.status,r.certified,r.dw,r.gam,r.ndv,r.m,r.t,r.msg);
if isappdata(0,'PF_ROW'),   rmappdata(0,'PF_ROW');  end
end
