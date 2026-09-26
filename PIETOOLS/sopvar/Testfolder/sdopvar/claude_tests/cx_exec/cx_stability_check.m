function R = cx_stability_check(cases)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = CX_STABILITY_CHECK(CASES) compares, WITHOUT SOLVING, the stock and the
% container program of each stability case (cx_cases_stability): the SDP
% shape of both (cx_shape, pre-solve; these LPIs have no 'ineq' constraints,
% so pre- and post-solve shapes coincide) and both assembly times. This is
% the whole comparison for the 2-D cases, which are never solved here.
%
% Also measures the two data-dependent quantities each transcription
% re-derives from containers, against the stock rule on the stock objects:
%   Q-form   Qdeg: cx_stability_lpivar_degs(P) vs get_lpivar_degs(Pop,Top);
%   2-D      eq_opts: cx_stability_eq_opts_2D(Q) vs get_eq_opts_2D(Qop),
%            the stock objects rebuilt by the stock executive's own lines.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(cases),   cases = [cases{:}];     end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
R = cell(1,numel(cases));
for ic = 1:numel(cases)
    c = cases(ic);
    st = cx_settings(c.setname,c.solver);
    PIE = cx_plant(c.plant{:});
    r = struct('id',c.id,'exec',c.exec,'err','','extra','');
    try
        t = tic;    prog0 = cx_stock_capture(c.exec,PIE,st);    r.tstock = toc(t);
        t = tic;    prog = feval(['cx_' c.exec],PIE,st);        r.tcx = toc(t);
        r.stock = cx_shape(prog0);  r.cx = cx_shape(prog);
        r.same = isequal(r.stock,r.cx);
        switch c.exec
            case {'PIE2PDEstability','PIE2PDEstability_dual'}
                [qs,qc] = qdeg_pair(PIE,st,strcmp(c.exec,'PIE2PDEstability_dual'));
                r.extra = sprintf('Qdeg stock %s cx %s',mat2str(full(qs)),mat2str(qc));
            case {'stability_2D','stability_dual_2D'}
                [es,ec] = eqopts_pair(PIE,st,strcmp(c.exec,'stability_dual_2D'));
                r.extra = sprintf('eq exclude stock %s cx %s | sep stock %s cx %s', ...
                    mat2str(double(es.exclude)),mat2str(double(ec.exclude)), ...
                    mat2str(double(es.sep)),mat2str(double(ec.sep)));
        end
    catch ME
        r.err = sprintf('[%s] %s (%s:%d)',ME.identifier,ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    if isempty(r.err)
        fprintf(['  %-14s %-22s same %d | stock ndv %d Kf %d m %d Ks %s (%.1fs) | '...
                 'cx ndv %d Kf %d m %d Ks %s (%.1fs) %s\n'],r.id,r.exec,r.same, ...
            r.stock.ndv,r.stock.Kf,r.stock.m,mat2str(r.stock.Ks),r.tstock, ...
            r.cx.ndv,r.cx.Kf,r.cx.m,mat2str(r.cx.Ks),r.tcx,r.extra);
    else
        fprintf('  %-14s %-22s ERROR %s\n',r.id,r.exec,r.err);
    end
    R{ic} = r;
end
R = [R{:}];
end

function [qs,qc] = qdeg_pair(PIE,st,dual)
% Stock lines 112-123 (dual: 115-128) on opvars, and the container P.
Top = PIE.T;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pop] = poslpivar(prog,Top.dim,st.dd1,st.options1);
if st.override1~=1,   [prog,P2] = poslpivar(prog,Top.dim,st.dd12,st.options12);   Pop = Pop+P2;   end
Tm = opvar2copvar(Top);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
if dual,    side = 'out';   else,   side = 'in';    end
[prog,Pm] = cx_stability_pos(prog,Tm,side,PIE.dom,st.dd1,st.options1);
if st.override1~=1,   [~,P2m] = cx_stability_pos(prog,Tm,side,PIE.dom,st.dd12,st.options12);   Pm = Pm+P2m;   end
if dual
    qs = get_lpivar_degs(Pop + st.eppos2*Top*Top',Top');
    qc = cx_stability_lpivar_degs(Pm + st.eppos2*Tm*Tm');
else
    qs = get_lpivar_degs(Pop + st.eppos2*Top'*Top,Top);
    qc = cx_stability_lpivar_degs(Pm + st.eppos2*Tm'*Tm);
end
end

function [es,ec] = eqopts_pair(PIE,st,dual)
% PIETOOLS_stability_2D.m:137-191 (dual: 139-193), LF psatz terms off as in
% the shipped files, on opvar2d objects; and the container Q.
s2 = st.settings_2d;    Top = PIE.T;    Aop = PIE.A;    e = s2.eppos;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[~,Pop] = poslpivar_2d(prog,Top.dim,s2.LF_deg,s2.LF_opts);
np = Pop.dim(:,1);
Pop = Pop + opvar2d(blkdiag(e(1)*eye(np(1)),e(2)*eye(np(2)),e(3)*eye(np(3)),e(4)*eye(np(4))), ...
                    Pop.dim,PIE.dom,PIE.vars);
if dual,    X = (Top*Pop)*Aop';     else,   X = Aop'*(Pop*Top);     end
es = get_eq_opts_2D(clean_opvar(X'+X,1e-12),s2.eq_opts,1e-12);
Tm = opvar2d2copvar(Top);   Am = opvar2d2copvar(Aop);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
if dual,    side = 'in';    else,   side = 'out';   end
[~,Pm] = cx_stability_pos2d(prog,Tm,side,PIE.dom,s2.LF_deg,s2.LF_opts);
Pm = Pm + opvar2d2copvar(opvar2d(blkdiag(e(1)*eye(np(1)),e(2)*eye(np(2)),e(3)*eye(np(3)), ...
                         e(4)*eye(np(4))),[np,np],PIE.dom,PIE.vars));
if dual,    Xm = (Tm*Pm)*Am';       else,   Xm = Am'*(Pm*Tm);       end
ec = cx_stability_eq_opts_2D(Xm'+Xm,s2.eq_opts,1e-12);
end
