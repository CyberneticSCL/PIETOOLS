function pc_post(part)
% PC_POST(PART) checks a posteriori the solves that the harness gate
% rejected for numerical reasons (maintainer 11/21/2026: "fix the storage
% output of the previous result and just find the slack"; the projection
% check: "you can try it").
%   Stage 1     the row as the harness ran it - the executive with the
%               tier, mods and form of the row - keeping the SDP data of
%               its last pass (post/sossolve.m, appdata SOS_SDP) and the
%               internals of the executive (info.aux, the instrument of
%               this worktree).
%   Stage 2     the storage fixed at its stage-1 value and only the slack
%               sought, sized as the stage-1 slack (its degree
%               specification and terms). 'P', 'Pd': every Gram block of
%               the stage-1 point clipped to the PSD cone, so that
%               P >= eppos I holds exactly. 'Q': the free part of the
%               stage-1 point (Q) projected onto the rows that involve Q
%               alone, T'*Q = Q'*T, which the derivative formula needs
%               exactly; then the slacks of P - (eppos2/2)*T'*T >= 0 (half
%               the stage-1 margin) and -D >= 0.
%   Projection  x_hat = x - At*((At'*At)\(At'*x - b)) on the SDP data of a
%               stage; lambda_min of every Gram block of x_hat. All of them
%               >= 0: x_hat solves the equalities and lies in the cone, up
%               to rounding.
% A stage passes the harness gate when the solver exit is 'certified'
% (|fr-1| <= 0.3, numerr 0, no pinf, no dinf), the executive status is
% 'optimal' and eta <= 1e-7. Soundness rows keep the stage-1 storage and
% seek the slack of the derivative of the PIE at an unstable parameter.
%
% INPUT
%   part:   'ctl' (controls), 'a', 'a2', 'b' (stability rows), 'c'
%           (library case 50, dual H2 'c'), 's' (SeDuMi rows); the rows
%           below.
% results_POST_<PART>.csv in this folder; no row starts after 18 minutes.
% MMP, 11/21/2026.
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
addpath(fullfile(here,'gate'));
addpath(fullfile(here,'post'));         % sossolve that keeps the SDP data
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY','SOS_SDP'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
% One resolution of the two instrumented files
w = which('sossolve');
if ~contains(w,[filesep,'post',filesep]),  error('pc_post:path','sossolve resolves to %s.',w);  end
w = which('PIETOOLS_stability_sop');
if ~contains(fileread(w),'wt1121v'),  error('pc_post:path','%s lacks the instrument.',w);  end
fprintf('## PC_POST %s %s\n## sossolve: %s\n## executive: %s\n',part,char(datetime('now')),which('sossolve'),w);
tstop = datetime('now') + minutes(18);
of = fullfile(here,['results_POST_',part,'.csv']);
fid = fopen(of,'w');
fprintf(fid,['id,system,param,form,tier,mods,s1_status,s1_fr,s1_exst,s1_eta,s1_ndec,s1_m,s1_lmin,', ...
             'p1_eta,p1_lmin,p1_lrel,p1_dx,prep,s2_status,s2_fr,s2_exst,s2_eta,s2_ndec,s2_m,s2_lmin,', ...
             'p2_eta,p2_lmin,p2_lrel,verdict,sec,note\n']);
fclose(fid);
SP = 'PIETOOLS_stability_sop';
% Rows {id, system, parameter, executive, form, tier, mods, unstable parameter}
switch part
    case 'ctl'
        J = {'C41','41x',[],SP,'Q','light','st.eppos = 1e-2',[];
             'C44','44min',0.4,SP,'Pd','stripped','',0.8};
    case 'a'
        J = {'R41','41x',[],SP,'Q','light','',[];
             'R44a','44min',0.4,SP,'Pd','light','st.epneg = 1e-4',[];
             'R44b','44min',0.4,SP,'Pd','stripped','st.epneg = 1e-3',[];
             'RX3a','x3min',0.375,SP,'Pd','heavy','st.epneg = 1e-3',[];
             'RX3b','x3min',0.4,SP,'Pd','heavy','st.epneg = 1e-3',[];
             'R36','PDE:36',[],SP,'Q','light','st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1',[]};
    case 'a2'       % the 'Pd' rows of 'a' again (projection to convergence), and 'Qd'
        J = {'R44a','44min',0.4,SP,'Pd','light','st.epneg = 1e-4',[];
             'R44b','44min',0.4,SP,'Pd','stripped','st.epneg = 1e-3',[];
             'RX3a','x3min',0.375,SP,'Pd','heavy','st.epneg = 1e-3',[];
             'RX3b','x3min',0.4,SP,'Pd','heavy','st.epneg = 1e-3',[];
             'R44q','44min',0.4,SP,'Qd','light','',[]};
    case 'b'
        J = {'R69','69mix2',[],SP,'Q','light','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]',[];
             'R36q','PDE:36',[],SP,'Q','light','st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1; st.epneg = 1e-3; st.sop.deg_Q{2} = [0 1 1]',[]};
    case 's'        % the SeDuMi rows of DIAGSED (pc_diagsed.m): both stages with SeDuMi 1.3
        addpath(genpath('C:/YALMIP/yalmip_install/tbxmanager/toolboxes/sedumi/1.3/pcwin64/sedumi_1_3_pcwin64'));
        fprintf('## sedumi: %s\n',which('sedumi'));
        sd = 'st.sos_opts.solver = ''sedumi''';
        J = {'S44','44min',0.4,SP,'Pd','light',sd,[];
             'S44e','44min',0.4,SP,'Pd','light',[sd,'; st.epneg = 1e-1'],[];
             'SX3','x3min',0.4,SP,'Pd','light',sd,[];
             'SX3e','x3min',0.4,SP,'Pd','light',[sd,'; st.epneg = 1e-2'],[]};
    case 'c'        % library case 50, the dual H2 test 'c': the stripped and light starts of EXS_50d2
        H2 = 'PIETOOLS_H2_norm_sop';
        J = {'H50a','PDE:50/H2_norm_dual',[],H2,'c','stripped','st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''',[];
             'H50b','PDE:50/H2_norm_dual',[],H2,'c','stripped','st.lpi.ineq.dD = 0',[];
             'H50c','PDE:50/H2_norm_dual',[],H2,'c','light','st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''',[]};
    otherwise
        error('pc_post:part','Unknown part %s.',part);
end
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
for i = 1:size(J,1)
    [id,sys,p,fn,form,tier,mods,pu] = J{i,:};
    if datetime('now') >= tstop
        fprintf('\n## PC_POST stop time %s reached\n',char(tstop));  return
    end
    fprintf('\n##### POST %s: %s %s %s %s [%s]  %s\n',id,sys,num2str(p),form,tier,mods,char(datetime('now')));
    t0 = tic;   note = '';  verdict = 'error';  pre = '';
    r1 = struct();  p1 = struct();  r2 = struct();  p2 = struct();
    info = [];  prog = [];  sdp1 = [];  S = [];  PIE = [];
    try
        PIE = build(sys,p,Cat);
        slv = 'mosek';  if contains(mods,'sedumi'),  slv = 'sedumi';  end   % as pc_diagsed.m
        S = lpisettings_sop(tier,0,'',slv);
        if ~isempty(mods),  eval([strrep(mods,'st.','S.'),';']);  end
        [r1,prog,info,sdp1] = stage1(PIE,S,fn,form);
        p1 = projcheck(sdp1);
        [r2,sdp2,pre] = stage2(PIE,S,form,prog,info,sdp1,[]);
        p2 = projcheck(sdp2);
        verdict = verdictof(r2,p2);
    catch ME
        note = oneline(ME);
        fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off'));
    end
    wrow(of,id,sys,p,form,tier,mods,r1,p1,pre,r2,p2,verdict,toc(t0),note);
    fprintf('##### POST %s: %s\n',id,verdict);
    % Soundness: the stage-1 storage against the PIE at the unstable parameter
    if isempty(pu) || isempty(info),  continue,  end
    t0 = tic;   note = '';  verdict = 'error';  pre = '';  r2 = struct();  p2 = struct();
    try
        PIEu = build(sys,pu,Cat);
        [r2,sdp2,pre] = stage2(PIE,S,form,prog,info,sdp1,PIEu);
        p2 = projcheck(sdp2);
        v = verdictof(r2,p2);
        if startsWith(v,'verified'),  verdict = ['UNSOUND - ',v];
        else,                         verdict = ['fails as required - ',v];
        end
    catch ME
        note = oneline(ME);
        fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off'));
    end
    wrow(of,[id,'u'],sys,pu,form,tier,mods,r1,p1,pre,r2,p2,verdict,toc(t0),note);
    fprintf('##### POST %su: %s\n',id,verdict);
end
fprintf('\n## PC_POST %s DONE %s\n',part,char(datetime('now')));
end


% ========================================================================
function [r,prog,info,sdp] = stage1(PIE,S,fn,form)
% The row as the harness ran it, with every output of the executive and
% the SDP data of its last pass.
if isappdata(0,'SOS_SDP'),  rmappdata(0,'SOS_SDP');  end
setappdata(0,'SOS_KEEP_SDP',true);
outs = cell(1,nargout(fn));
[outs{:}] = feval(fn,PIE,S,form);
prog = outs{1};     info = outs{end};
sdp = [];   if isappdata(0,'SOS_SDP'),  sdp = getappdata(0,'SOS_SDP');  end
r = solrec(prog);
r.exst = info.status;   r.eta = NaN;
if isfield(info,'hist') && ~isempty(info.hist),  r.eta = double(info.hist(end).eta);  end
fprintf('##### stage 1: %s fr %.4g exst %s eta %.3g ndec %d m %d\n',r.status,r.fr,r.exst,r.eta,r.ndec,r.m);
end

function [r,sdp,pre] = stage2(PIE,S,form,prog,info,sdp1,PIEu)
% The storage of stage 1 fixed, the slack sought. PIEu nonempty: the
% derivative of PIEu with that storage ('P', 'Pd').
if isempty(sdp1),  error('pc_post:nosdp','No SDP data from stage 1.');  end
if isfield(info,'gam'),  [r,sdp,pre] = stage2h2(PIE,S,form,prog,info,sdp1);  return,  end
if ~any(strcmp(form,{'P','Pd','FP','Q','Qd'})),  error('pc_post:form','No stage 2 for form %s.',form);  end
aux = info.aux;
% The storage (P; Q in 'Q') projected onto the rows of stage 1 that touch
% it alone: coefficients of D no slack monomial reaches, and T'*Q = Q'*T.
% A fixed storage off those rows leaves D with components outside the span
% of the slack, which the slack program cannot match (sossolve declares
% such a row infeasible before the solve). Gram blocks of P alternate with
% clipping to the PSD cone.
isS = storidx(prog,sdp1,{aux.Pm,aux.Dm});     % P, and Q in 'Q', 'Qd'
[x,sp] = projstor(sdp1,isS);
pre = sprintf('storage rows %d; eta %.2g to %.2g; inconsistency %.2g; storage lambda_min %.3g; clipped %.2g', ...
              sp.n,sp.eta0,sp.eta1,sp.ls,sp.lmin,sp.clip);
prog.solinfo.RRx = (sdp1.RR*x)*sdp1.bscl;
if isempty(PIEu)
    PIEc = exec_pie_sop(PIE);
    Dn = getsol_lpivar_sop(prog,aux.Dm);
    op = {optsli(aux.li)};
else
    PIEc = exec_pie_sop(PIEu);
    Pn = getsol_lpivar_sop(prog,aux.Pm);
    Tm = PIEc.T;    Am = PIEc.A;    ep = info.epneg;
    switch form
        case 'P'
            Dn = Tm'*Pn*Am + Am'*Pn*Tm;     if ep~=0,  Dn = Dn + 2*ep*(Tm'*Pn*Tm);  end
        case 'Pd'
            Dn = Tm*Pn*Am' + Am*Pn*Tm';     if ep~=0,  Dn = Dn + 2*ep*(Tm*Pn*Tm');  end
        otherwise
            error('pc_post:form','No soundness row for form %s.',form);
    end
    % The slack sized from the support of this derivative (lift rule, dw 2,
    % no pruning), so that its monomials are reachable: a failure is then
    % the solver's verdict on the inequality
    op = {struct('dw',2,'margin',0,'prune',false)};
    pre = [pre,'; derivative of the unstable PIE, slack sized from its support'];
end
tg = {-Dn};
if isempty(PIEu) && any(strcmp(form,{'P','Pd','Qd'}))
    % The margin of the stage-1 storage, P - eppos I >= 0 ('P', 'Pd') or
    % P - eppos2 T*T' >= 0 ('Qd'), by a slack sized from its own support:
    % the projection above can leave the Gram block of P indefinite
    st = exec_sop_settings(S);  tl = exec_tools_sop(PIEc,st);  st.sop = tl.sop;
    Pn = getsol_lpivar_sop(prog,aux.Pm);    ep = info.eppos;
    switch form
        case 'P',   tpos = Pn - tl.eye(1,ep);
        case 'Pd',  tpos = Pn - tl.eye(2,ep);
        case 'Qd',  tpos = Pn - ep(2)*(PIEc.T*PIEc.T');
    end
    tg = [{tpos}, tg];
    op = [{struct('dw',1,'margin',0,'prune',false)}, op];
    pre = sprintf('%s; P - margin >= 0 checked',pre);
end
if strcmp(form,'Q')
    % P = T'*Q >= eppos2 T'*T, the stage-1 margin; its slack sized from its
    % own support as in stage 1 (raise dw of the last pass, base 0), no
    % pruning
    st = exec_sop_settings(S);  tl = exec_tools_sop(PIEc,st);  st.sop = tl.sop;
    Pn = getsol_lpivar_sop(prog,aux.Pm);
    e2 = min(info.eppos(min(2,end):end));
    dws = info.dw;  if strcmp(st.sop.raise_Q,'deg'),  dws = 0;  end
    opos = struct('dw',dws,'margin',st.sop.margin,'cap',tl.capfor(Pn),'prune',false);
    tg = [{Pn - e2*(PIEc.T'*PIEc.T)}, tg];
    op = [{opos}, op];
    pre = sprintf('%s; P >= %.3g T''T',pre,e2);
end
[r,sdp,li] = slackonly(PIEc,S,tg,op);
pre = sprintf('%s; slack entries %s (stage 1: %d)',pre,mat2str(cellfun(@(s) s.ndec,li)),aux.li.ndec);
end

function [r,sdp,pre] = stage2h2(PIE,S,form,prog,info,sdp1)
% PIETOOLS_H2_norm_sop, form 'c': R, Q, W and gam fixed at their stage-1
% values (R, W Gram blocks clipped, all projected onto the rows that touch
% them alone, among them T*Q = R), and the slacks of
%   [-g I, B1'; B1, Q'*A' + A*Q] <= 0,   [W, C1*Q; Q'*C1', R] >= 0
% sought, g = max(gam, trace W) (both constraints hold for any g at least
% as large); sized as stage 1, no pruning.
if ~strcmp(form,'c'),  error('pc_post:form','No H2 stage 2 for form %s.',form);  end
aux = info.aux;     gv = aux.obj;
isS = storidx(prog,sdp1,{aux.Rm,aux.Qm,aux.Wm,gv});
[x,sp] = projstor(sdp1,isS);
prog.solinfo.RRx = (sdp1.RR*x)*sdp1.bscl;
Rn = getsol_lpivar_sop(prog,aux.Rm);    Qn = getsol_lpivar_sop(prog,aux.Qm);    Wn = getsol_lpivar_sop(prog,aux.Wm);
dtab = prog.decvartable;    if isstring(dtab),  dtab = cellstr(dtab);  end
gn = prog.solinfo.RRx(find(strcmp(dtab,gv.dvarname{1}),1));
PIEc = exec_pie_sop(PIE);
st = exec_sop_settings(S);  tl = exec_tools_sop(PIEc,st);  st.sop = tl.sop;
tr = trace_rn_sop(Wn);
if isa(tr,'dpvar'),  trn = full(tr.C(1));  else,  trn = double(tr);  end
gc = max(gn,trn);
Tm = PIEc.T;    Am = PIEc.A;    Bw = PIEc.B1;   Cz = PIEc.C1;   Iw = tl.ident(PIEc.B1,2); %#ok<NASGU>
Dneg = [-(gc*Iw), Bw'; Bw, Qn'*Am' + Am*Qn];
Dpos = [Wn, Cz*Qn; Qn'*Cz', Rn];
opn = optsli(aux.li);
opp = struct('like',Rn,'dw',st.sop.dw_Q+info.dw,'margin',st.sop.margin,'cap',tl.capfor(Dpos),'prune',false);
[r,sdp,li] = slackonly(PIEc,S,{-Dneg, Dpos},{opn, opp});
pre = sprintf(['storage rows %d; eta %.2g to %.2g; inconsistency %.2g; Gram lambda_min %.3g; clipped %.2g; ', ...
               'gam %.8g; trace W %.8g; slack entries %s (stage 1 negativity: %d)'], ...
              sp.n,sp.eta0,sp.eta1,sp.ls,sp.lmin,sp.clip,gn,trn,mat2str(cellfun(@(s) s.ndec,li)),aux.li.ndec);
end

function [r,sdp,li] = slackonly(PIEc,S,tg,op)
% A program of slacks only: tg{k} >= 0 with the options op{k}, solved.
st = exec_sop_settings(S);  tl = exec_tools_sop(PIEc,st);  st.sop = tl.sop;
prog = lpiprogram(PIEc.vars(:,1),PIEc.vars(:,2),PIEc.dom,st);
li = cell(size(tg));
for k = 1:numel(tg)
    [prog,~,li{k}] = lpi_ineq_sop(prog,tg{k},op{k});
end
if isappdata(0,'SOS_SDP'),  rmappdata(0,'SOS_SDP');  end
setappdata(0,'SOS_KEEP_SDP',true);
prog = lpisolve(prog,st.sos_opts);
cl = lpi_classify_sop(prog,[]);
sdp = [];   if isappdata(0,'SOS_SDP'),  sdp = getappdata(0,'SOS_SDP');  end
r = solrec(prog);   r.exst = cl.status;     r.eta = cl.eta;
fprintf('##### stage 2: %s fr %.4g exst %s eta %.3g ndec %d m %d\n',r.status,r.fr,r.exst,r.eta,r.ndec,r.m);
end

function o = optsli(li)
% The options that rebuild the slack of LI: its degree specification
% (with the cap subsets) and its terms.
% No pruning: pruning reads the nonzero cells of the numeric target, which
% kept fewer basis operators than stage 1 (1,352 against 1,568 entries on
% the 41 control); unpruned, the slack contains the stage-1 slack.
o = struct('deg',{li.deg},'psatz',li.terms.codes,'psatz_offset',li.terms.offsets,'margin',li.margin,'prune',false);
if isfield(li.terms,'rows') && ~isempty(li.terms.rows)
    o.psatz_rows = li.terms.rows;   o.psatz_rows_offset = li.terms.rows_offsets;
end
end

% ========================================================================
function r = solrec(prog)
% The solver exit and the verdict of LB_EXEC.
r = struct('status','nosolinfo','fr',NaN,'pinf',NaN,'dinf',NaN,'numerr',NaN,'ndec',NaN,'m',NaN);
if isfield(prog,'decvartable'),  r.ndec = numel(prog.decvartable);  end
if isfield(prog,'solinfo') && isfield(prog.solinfo,'info')
    in = prog.solinfo.info;
    if isfield(in,'feasratio'),  r.fr = double(in.feasratio);  end
    for k = {'pinf','dinf','numerr'}
        if isfield(in,k{1}),  r.(k{1}) = double(in.(k{1}));  end
    end
    if isfield(prog.solinfo,'y'),  r.m = numel(prog.solinfo.y);  end
end
if any(isnan([r.fr r.pinf r.dinf r.numerr])),                    r.status = 'nosolinfo';
elseif abs(r.fr-1)<=0.3 && ~r.numerr && ~r.pinf && ~r.dinf,      r.status = 'certified';
elseif abs(r.fr-1)<=0.3 && r.numerr && ~r.pinf && ~r.dinf,       r.status = 'inaccurate';
elseif r.pinf || r.dinf || abs(r.fr+1)<=0.1,                     r.status = 'infeasible';
else,                                                            r.status = 'unknown';
end
end

function v = verdictof(r2,p2)
gate  = strcmp(f(r2,'status'),'certified') && strcmp(f(r2,'exst'),'optimal') && g(r2,'eta') <= 1e-7;
exact = ~isnan(g(p2,'lmin')) && g(p2,'lmin') >= 0;
if exact && gate,   v = 'verified';
elseif exact,       v = 'verified (projection only)';
elseif gate,        v = 'stage-2 gate only';
else,               v = 'not verified';
end
end

% ========================================================================
function p = projcheck(sdp)
% lambda_min of the Gram blocks before and after the projection of x onto
% the affine set At'*x = b.
p = struct('lmin0',NaN,'lmin',NaN,'lrel',NaN,'eta0',NaN,'eta1',NaN,'dx',NaN);
if isempty(sdp),  return,  end
At = sdp.At;    b = full(sdp.b);    x = full(sdp.x);    K = sdp.K;
n = getk(K,'f') + getk(K,'l') + sum(K.s.^2);
if n~=numel(x),  error('pc_post:K','K describes %d entries, x has %d.',n,numel(x));  end
p.lmin0 = min(gramlmin(x,K));
p.eta0 = etaof(At,b,x);
xh = projaff(At,b,x);
p.eta1 = etaof(At,b,xh);
[lam,lmax] = gramlmin(xh,K);
p.lmin = min(lam);
ok = lmax>0 & isfinite(lam);
p.lrel = min(lam(ok)./lmax(ok));
p.dx = norm(xh-x)/max(norm(x),realmin);
fprintf('##### projection: eta %.3g to %.3g; lambda_min %.3g to %.3g (relative %.3g); |dx|/|x| %.3g\n', ...
        p.eta0,p.eta1,p.lmin0,p.lmin,p.lrel,p.dx);
end

function xh = projaff(At,b,x)
% The point of At'*x = b nearest to x. Up to 5e7 entries: minimum-norm
% least-squares steps on the dense At' (complete orthogonal decomposition,
% LSQMINNORM), which resolve a rank-deficient, ill-conditioned row set to
% rounding; there the normal equations stalled at eta 5e-9 (44min 'Pd').
% Larger: x - At*((At'*At)\r) with refinement.
if numel(At)<=5e7
    xh = x;     A = full(At');
    for it = 1:3
        xh = xh - lsqminnorm(A,full(A*xh - b));
    end
    return
end
M = At'*At;     m = size(M,1);
[R,flag,pp] = chol(M,'vector');             % M(pp,pp) = R'*R
if flag~=0
    mu = 1e-12*max(full(diag(M)));          % rank-deficient row set
    [R,flag,pp] = chol(M+mu*speye(m),'vector');
    if flag~=0,  error('pc_post:chol','At''*At does not factor (%d rows).',m);  end
end
% Refinement until the backward error stops falling (30 steps at most): the
% regularized factor converges slowly along small singular values
xh = x;     e0 = inf;
for it = 1:30
    r = full(At'*xh - b);
    w = zeros(m,1);     w(pp) = R\(R'\r(pp));
    xh = xh - At*w;
    e1 = etaof(At,b,xh);
    if e1<1e-15 || e1>0.5*e0,  break,  end
    e0 = e1;
end
end

function c = lsres(At,b)
% Inconsistency of At'*z = b: the part of b along the eigenvectors of
% At'*At with eigenvalues below 1e-12 of the largest, relative to |b|
% (0: consistent). Dense; NaN above 3000 rows.
c = NaN;
if size(At,2)>3000,  return,  end
[V,E] = eig(full(At'*At));  e = diag(E);
sm = e < 1e-12*max(e);
c = norm(V(:,sm)'*b)/max(norm(b),realmin);
end

function e = etaof(At,b,x)
% The row-normwise backward error of SOSSOLVE.
r = full(At'*x - b);
rn = full(sqrt(sum(At.^2,1)))';     bb = abs(full(b));
on = rn>1e-12*max(rn) | bb>1e-12*max(bb);
e = max([0; abs(r(on))./(rn(on)*max(abs(x))+bb(on)+realmin)]);
end

function [lam,lmax] = gramlmin(x,K)
% Smallest and largest |eigenvalue| of every Gram block (last entry: the
% linear part).
off = getk(K,'f') + getk(K,'l');
ns = numel(K.s);    lam = inf(ns+1,1);  lmax = zeros(ns+1,1);
for i = 1:ns
    s = K.s(i);     idx = off+(1:s*s);  off = off + s*s;
    if s==0,  continue,  end
    X = reshape(full(x(idx)),s,s);  X = (X+X')/2;
    e = eig(X);     lam(i) = min(e);    lmax(i) = max(abs(e));
end
nl = getk(K,'l');
if nl>0
    v = x(getk(K,'f')+(1:nl));  lam(end) = min(v);  lmax(end) = max(abs(v));
end
end

function [xc,cl] = clipgram(x,K)
% Every Gram block of x projected onto the PSD cone (negative eigenvalues
% set to 0); cl the largest |eigenvalue| removed.
xc = full(x);   cl = 0;
nf = getk(K,'f');   nl = getk(K,'l');   off = nf + nl;
if nl>0
    v = xc(nf+(1:nl));  cl = max(cl,-min([v;0]));   xc(nf+(1:nl)) = max(v,0);
end
for i = 1:numel(K.s)
    s = K.s(i);     idx = off+(1:s*s);  off = off + s*s;
    if s==0,  continue,  end
    X = reshape(xc(idx),s,s);   X = (X+X')/2;
    [V,E] = eig(X);     e = diag(E);
    if min(e)<0
        cl = max(cl,-min(e));
        X = V*diag(max(e,0))*V';    X = (X+X')/2;
    end
    xc(idx) = X(:);
end
end

function isS = storidx(prog,sdp,ops)
% The entries of the SDP point x that feed the decision variables of the
% operators OPS (a cell of 'cdopvar' or 'dpvar' objects), through
% RRx = RR*x with rows in decvartable order.
if ~iscell(ops),  ops = {ops};  end
Zd = {};
for k = 1:numel(ops)
    if isempty(ops{k}),  continue,  end
    if isa(ops{k},'dpvar'),  z = ops{k}.dvarname;  else,  z = ops{k}.Zd;  end
    if isstring(z),  z = cellstr(z);  end
    Zd = [Zd; z(:)]; %#ok<AGROW>
end
dtab = prog.decvartable;    if isstring(dtab),  dtab = cellstr(dtab);  end
rows = find(ismember(dtab(:),Zd(:)));
isS = full(any(sdp.RR(rows,:)~=0,1))';
if ~any(isS),  error('pc_post:storage','No SDP entry feeds the storage.');  end
end

function [x,sp] = projstor(sdp,isS)
% The storage entries of x projected onto the rows that touch them alone,
% alternating with the PSD projection of the storage Gram blocks; the last
% step is affine, so those rows hold up to rounding.
At = sdp.At;    b = full(sdp.b);    x = full(sdp.x);    K = sdp.K;
tS = full(sum(At(isS,:)~=0,1))'>0;
tO = full(sum(At(~isS,:)~=0,1))'>0;
rows = find(tS & ~tO);
idx = find(isS);
% Gram blocks lying in the storage
off = getk(K,'f') + getk(K,'l');    B = zeros(0,2);
for i = 1:numel(K.s)
    s = K.s(i);
    if s>0 && all(isS(off+(1:s*s))),  B(end+1,:) = [off s];  end %#ok<AGROW>
    off = off + s*s;
end
sp = struct('n',numel(rows),'eta0',NaN,'eta1',NaN,'lmin',NaN,'clip',0,'ls',NaN);
if ~isempty(rows)
    As = At(idx,rows);  bs = b(rows);
    sp.eta0 = etaof(As,bs,x(idx));
    sp.ls = lsres(As,bs);
    for it = 1:5
        x(idx) = projaff(As,bs,x(idx));
        cl = 0;
        for j = 1:size(B,1)
            ii = B(j,1)+(1:B(j,2)^2);
            X = reshape(x(ii),B(j,2),B(j,2));  X = (X+X')/2;
            [V,E] = eig(X);     e = diag(E);
            if min(e)<0,  cl = max(cl,-min(e));  X = V*diag(max(e,0))*V';  X = (X+X')/2;  end
            x(ii) = X(:);
        end
        sp.clip = max(sp.clip,cl);
        if cl==0,  break,  end
    end
    x(idx) = projaff(As,bs,x(idx));
    sp.eta1 = etaof(As,bs,x(idx));
end
lm = inf;
for j = 1:size(B,1)
    X = reshape(x(B(j,1)+(1:B(j,2)^2)),B(j,2),B(j,2));
    lm = min(lm,min(eig((X+X')/2)));
end
sp.lmin = lm;
fprintf('##### storage: %d rows touch it alone, eta %.3g to %.3g, inconsistency %.3g; %d Gram blocks, lambda_min %.3g, clipped %.2g\n', ...
        sp.n,sp.eta0,sp.eta1,sp.ls,size(B,1),sp.lmin,sp.clip);
end

function [x,sa] = projfree(sdp)
% The free part of x projected onto the rows that involve free variables
% only ('Q': T'*Q - Q'*T = 0).
At = sdp.At;    b = full(sdp.b);    x = full(sdp.x);    nf = getk(sdp.K,'f');
if nf==0,  error('pc_post:nofree','No free variables in the SDP.');  end
inG = full(sum(At(nf+1:end,:)~=0,1))'>0;
inF = full(sum(At(1:nf,:)~=0,1))'>0;
rows = find(inF & ~inG);
sa = struct('n',numel(rows),'eta0',NaN,'eta1',NaN);
if isempty(rows),  return,  end
Af = At(1:nf,rows);     bf = b(rows);   xf = x(1:nf);
sa.eta0 = etaof(Af,bf,xf);
xf = projaff(Af,bf,xf);
sa.eta1 = etaof(Af,bf,xf);
x(1:nf) = xf;
end

function v = getk(K,fld)
v = 0;  if isfield(K,fld) && ~isempty(K.(fld)),  v = sum(K.(fld));  end
end

% ========================================================================
function PIE = build(sys,p,Cat)
switch sys
    case '41x',                 PIE = build41(1,1);
    case {'44min','x3min'},     PIE = build44e(sys,p);
    case '69mix2',              PIE = build69m();
    otherwise                   % a catalog key
        c = Cat(find(strcmp({Cat.key},sys),1));
        if isempty(c),  error('pc_post:key','No catalog row %s.',sys);  end
        PIE = lb_build_system(c,p,true);
end
end

function PIE = build41(lam,tau)
% PDE 41 (PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input) with its delay
% as the transport x2 on s2 in [1, 1+tau] (pc_diag41b.m, 11/20/2026).
pvar s1 s2 t
clear stateNameGenerator
pde_var x1 x2 control u
x1.vars = s1;   x1.dom = [0,1];
x2.vars = s2;   x2.dom = [1,1+tau];
PDE_t = [diff(x1,t)==diff(x1,s1,2)+lam*x1;
         diff(x2,t)==diff(x2,s2);
         subs(x1,s1,0)==0;  subs(x1,s1,1)==subs(x2,s2,1);
         subs(x2,s2,1+tau)==u];
PIE = convert(PDE_t,'pie');
end

function PIE = build44e(v,mu)
% The minimal position form of PDE 44 (PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay):
% u (order 2), v = u_t (order 1), the delayed velocity z (order 1), on s
% ('44min') or with z on a second variable s2 ('x3min', the layout of
% PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D); pc_diag44e.m, 11/20/2026.
k = 1;  tau = 1;
pvar s s2 t
clear stateNameGenerator
u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);
if strcmp(v,'44min'),  z = pde_var(s,[0,1],1);   zs = s;
else,                  z = pde_var(s2,[0,1],1);  zs = s2;
end
PDE_t = [diff(u,'t')==v1;
         diff(v1,'t')==diff(u,s,2);
         diff(z,'t')==-(1/tau)*diff(z,zs);
         subs(u,s,0)==0;  subs(v1,s,0)==0;
         subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,zs,1);
         subs(z,zs,0)==subs(v1,s,1)];
PIE = convert(PDE_t,'pie');
end

function PIE = build69m()
% '69mix2', the mixed position-strain form of PDE 69
% (PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty); pc_diag69m.m, 11/20/2026.
pvar s1 s2 t
clear stateNameGenerator
u = pde_var([s1;s2],[0,1;0,1],[0;2]);
e1 = pde_var([s1;s2],[0,1;0,1],[1;2]);
vv = pde_var([s1;s2],[0,1;0,1],[1;2]);
w = pde_var([s1;s2],[0,1;0,1],[2;2]);
PDE_t = [diff(u,'t')==vv;
         diff(vv,'t')==diff(e1,s1)+diff(u,s2,2);
         diff(e1,'t')==diff(vv,s1);
         diff(w,'t')==diff(w,s1,2)+diff(w,s2,2);
         subs(u,s2,0)==0;           subs(u,s2,1)==0;
         subs(e1,s1,0)==-subs(diff(w,s1),s1,0);
         subs(e1,s2,0)==0;          subs(e1,s2,1)==0;
         subs(vv,s1,1)==0;          subs(vv,s2,0)==0;    subs(vv,s2,1)==0;
         subs(w,s1,1)==0;           subs(w,s1,0)==subs(vv,s1,0);
         subs(w,s2,0)==0;           subs(w,s2,1)==0];
PIE = convert(PDE_t,'pie');
end

% ========================================================================
function wrow(of,id,sys,p,form,tier,mods,r1,p1,pre,r2,p2,verdict,sec,note)
c = {id, sys, num2str(p), form, tier, ['"',strrep(mods,'''',''),'"'], ...
     f(r1,'status'), n4(g(r1,'fr')), f(r1,'exst'), n3(g(r1,'eta')), n0(g(r1,'ndec')), n0(g(r1,'m')), n3(g(p1,'lmin0')), ...
     n3(g(p1,'eta1')), n3(g(p1,'lmin')), n3(g(p1,'lrel')), n3(g(p1,'dx')), ['"',strrep(pre,'"',''''),'"'], ...
     f(r2,'status'), n4(g(r2,'fr')), f(r2,'exst'), n3(g(r2,'eta')), n0(g(r2,'ndec')), n0(g(r2,'m')), n3(g(p2,'lmin0')), ...
     n3(g(p2,'eta1')), n3(g(p2,'lmin')), n3(g(p2,'lrel')), verdict, sprintf('%.1f',sec), ['"',note,'"']};
fid = fopen(of,'a');  fprintf(fid,'%s\n',strjoin(c,','));  fclose(fid);
end

function s = oneline(ME)
s = ME.message;
if ~isempty(ME.stack),  s = sprintf('%s @ %s:%d',s,ME.stack(1).name,ME.stack(1).line);  end
s = strrep(strrep(strrep(s,',',';'),newline,' '),'"','''');
end

function s = n3(x),  s = sprintf('%.3g',x);  end
function s = n4(x),  s = sprintf('%.4g',x);  end
function s = n0(x),  s = sprintf('%.0f',x);  end
function s = f(r,k),  s = '';   if isstruct(r) && isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isstruct(r) && isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end
