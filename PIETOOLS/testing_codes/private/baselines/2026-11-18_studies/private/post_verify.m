function [ok,r2,pre] = post_verify(PIE,S,form,prog,info,delta)
% [OK,R2,PRE] = POST_VERIFY(PIE,S,FORM,PROG,INFO,DELTA) checks a solve of
% PIETOOLS_stability or PIETOOLS_H2_norm a posteriori: the storage
% operators of the solve are fixed at their values (getsol_lpivar on the
% operators the executive returns in info.aux) and only the slacks are
% sought (maintainer, 11/21/2026: "fix the storage output of the previous
% result and just find the slack"; "All we need is the storage operator
% itself, which is obtained from getsol").
%
% The slacks, in a program of slack variables only:
%   every stability form: -D >= 0, D the derivative operator of the solve;
%   the direct forms also the margin of P = sym(T'*Q), not positive by
%   construction, with the slack of the solve (info.aux.lip, Npm):
%   'Q' P - eppos2 T'*T >= 0; 'Qd' in two or more variables
%   (the direct form on (T',A')) P - eppos2 T*T' >= 0, P = T*Q; 'FQ' in two
%   or more variables P >= 0. In the storage forms ('P', 'Pd', 'FP'; 'Qd',
%   'FQ' in one variable) P is the Gram parametrization of the solve,
%   positive by its construction;
%   H2 'c'         [-g I, B1'; B1, Q'*A' + A*Q] <= 0 and [W, C1*Q; Q'*C1', R] >= 0,
%                  R, Q, W fixed, g = max(gam, trace W).
% -D is sought first with the slack of the solve (its degree specification,
% terms and kept basis operators, on the spaces of the solve), then, if that
% does not pass, with a slack sized from the support of -D itself one degree
% higher, which reaches every monomial of the fixed operator. With the slack
% of the solve, the coefficients of -D that no slack variable reaches are
% residue of the solve and are set to zero; the check fails if one exceeds
% the equality residual of the solve or 1e-7 of the largest |D| (H2 'c': if
% any is nonzero, since the executive returns no slack operator to measure
% that residual).
% The forms with a storage equality ('Q', 'Qd', 'FQ': T'*Q = Q'*T or T*Q =
% P; H2 'c': T*Q = R) need it to hold for the derivative formula; the fixed
% operators satisfy it to the residual of the solve, so the check applies
% only when the solve has eta <= 1e-7 (else OK is false, PRE says why).
% With DELTA > 0 the derivative is that of the shifted PIE T x_t = (A +
% DELTA*T) x, whose spectrum is the stable one moved right by DELTA:
% D + 2*DELTA*W, W the storage seen by the state (T'PT, TPT', or P), a
% soundness row that a sound check must fail once DELTA exceeds the decay rate.
%
% INPUT
% - PIE:    'pie_struct' object, the PIE of the solve;
% - S:      'struct', the settings of the solve;
% - form:   char, the form of the solve;
% - prog:   'struct', the solved program (the first output of the executive);
% - info:   'struct', the last output of the executive, with the field aux;
% - delta:  (optional) double >= 0, the shift of the soundness row (0).
% OUTPUT
% - ok:     logical, a slack-only program passes the harness gate (solver
%           exit 'certified', executive status optimal, eta <= 1e-7);
% - r2:     'struct', its exit: status, fr, pinf, dinf, numerr, ndec, m,
%           exst, eta, and gam for H2;
% - pre:    char, how the check went (slack sizes, the attempt that passed).
% MMP, 11/21/2026.
if nargin<6 || isempty(delta),  delta = 0;  end
if ~isfield(info,'aux'),  error('post_verify:aux','The executive returned no info.aux.');  end
r2 = struct('status','not run','exst','','eta',NaN,'fr',NaN,'ndec',NaN,'m',NaN);
eta1 = NaN;
if isfield(info,'hist') && ~isempty(info.hist),  eta1 = double(info.hist(end).eta);  end
qtype = any(strcmp(form,{'Q','Qd','FQ','c'}));
if qtype && ~(eta1 <= 1e-7)
    ok = false;
    pre = sprintf('not checked: form %s needs its storage equality, and the solve has eta %.2g > 1e-7',form,eta1);
    return
end
if isfield(info,'gam') && strcmp(form,'c')
    [ok,r2,pre] = stage2h2(PIE,S,prog,info);
else
    [ok,r2,pre] = stage2(PIE,S,form,prog,info,delta);
end
end

% ------------------------------------------------------------------------
function [ok,r,pre] = stage2(PIE,S,form,prog,info,delta)
% The stability forms: the storage fixed (getsol), the slacks sought.
if ~any(strcmp(form,{'P','Pd','FP','Q','Qd','FQ'}))
    error('post_verify:form','No a posteriori check for form %s.',form);
end
aux = info.aux;     nv = info.nv;
direct = strcmp(form,'Q') || (nv>=2 && any(strcmp(form,{'Qd','FQ'})));
PIEc = exec_pie(PIE);
st = exec_sop_settings(S);  tl = exec_tools(PIEc,st);  st.sop = tl.sop;
Pn = getsol_lpivar(prog,aux.Pm);
Dn = getsol_lpivar(prog,aux.Dm);
Tm = PIEc.T;    if direct && strcmp(form,'Qd'),  Tm = Tm';  end     % 'Qd' direct: the pair (T',A')
ep = info.eppos;    e2 = min(ep(min(2,end):end));
pre = '';
if delta>0                          % the shifted PIE: D + 2*delta*W
    switch form
        case {'P','FP'},    W = PIEc.T'*Pn*PIEc.T;
        case 'Pd',          W = PIEc.T*Pn*PIEc.T';
        otherwise,          W = Pn;     % P = T'*Q (T*Q): (A+dT)'Q + Q'(A+dT) = D + 2d P
    end
    Dn = Dn + (2*delta)*W;
    pre = sprintf('shifted by %g; ',delta);
end
% The margin of the storage: in the direct forms P = sym(T'*Q) is not
% positive by construction, so the slack of its margin in the solve
% (info.aux.lip), rebuilt exactly, checks it, with its residue rows zeroed
% as for -D below (a slack sized from the numeric P differs from the
% solve's and failed case 41, whose -D passed); in the storage forms ('P',
% 'Pd', 'FP', and 'Qd', 'FQ' in 1-D) P is the Gram parametrization
% Z'*G*Z (+ eppos I) of the solve, positive by its construction, and only
% -D is sought
tg0 = {};   op0 = {};   zt0 = [];
if direct
    if strcmp(form,'FQ'),   tpos = Pn;      else,   tpos = Pn - e2*(Tm'*Tm);   end
    tg0 = {tpos};
    if isfield(aux,'lip') && ~isempty(aux.lip)
        op0 = {struct('like_li',aux.lip,'ref',aux.Pm)};
        Rp = tpos - getsol_lpivar(prog,aux.Npm);       % its residual in the solve
        if aux.lip.margin>0,  Rp = Rp + aux.lip.margin*aux.lip.eD;  end
        tmax = maxabs_opvar(tpos);
        zt0 = min(maxabs_opvar(Rp) + 1e-12*tmax, 1e-7*tmax);        % as ztol for -D below
        pre = sprintf('%sP margin: residual of the solve %.2g, eppos2 %.3g; ',pre,maxabs_opvar(Rp),e2);
    else                                % an executive without aux.lip: sized as the solve sized it
        dws = info.dw;  if strcmp(st.sop.raise_Q,'deg'),  dws = 0;  end
        op0 = {struct('dw',dws,'margin',st.sop.margin,'cap',tl.capfor(Pn))};
    end
end
% The coefficients of D that the slack of the solve cannot reach were driven
% to zero by equality rows in the storage variables alone, to the residual
% of the solve; fixed, they are rows "0 = residue" that make the slack
% program inconsistent. slack_like zeroes exactly those coefficients (rows
% of P - Pop == 0 without a slack variable) and keeps every other one, however
% small: a small reachable coefficient is a Gram diagonal entry supporting
% larger off-diagonal ones (zeroing by magnitude changed verdicts both ways,
% cases 36 and 62, 11/21). On those rows the residue equals the residual of
% the solve, R = -D + margin*eD - N (N its slack, getsol); a zeroed value
% above max|R| (plus round-off) is not residue, and the check fails, as it
% does above 1e-7 of the largest |D| (the tolerance of the gate): an
% inaccurate solve (case NDS11, max|R| 1.4e-4 of max|D|) would otherwise
% admit zeroing a part of D that matters
Nn = getsol_lpivar(prog,aux.Nm);
Rn = -getsol_lpivar(prog,aux.Dm) - Nn;
if isfield(aux.li,'margin') && aux.li.margin>0,  Rn = Rn + aux.li.margin*aux.li.eD;  end
rmax = maxabs_opvar(Rn);   dmax = maxabs_opvar(Dn);
ztol = min(rmax + 1e-12*dmax, 1e-7*dmax);
pre = sprintf('%sresidual of the solve %.2g, largest D %.2g; ',pre,rmax,dmax);
% -D: the slack of the solve, then one sized from the support of -D, one
% degree higher (reaches every monomial of the fixed operator)
opA = struct('like_li',aux.li,'ref',aux.Dm);     % the slack of the solve, rebuilt exactly
opB = struct('dw',info.dw+1,'margin',0);
[r,li] = slackonly(PIEc,S,[tg0,{-Dn}],[op0,{opA}]);
dz = li{end}.drop;     okz = dz<=ztol;
if ~isempty(zt0),  dz = [li{1}.drop, dz];  okz = okz && dz(1)<=zt0;  end   % [P margin, -D]
ok = gatepass(r) && okz;   att = sprintf('slack of the solve, residue %s zeroed',mat2str(dz,2));
rA = sprintf('first: %s fr %.3g eta %.2g, residue %s',r.status,r.fr,r.eta,mat2str(dz,2));
if ~okz,  rA = sprintf('%s above the residual of the solve: not residue',rA);  end
if ~ok
    try
        [rB,liB] = slackonly(PIEc,S,[tg0,{-Dn}],[op0,{opB}]);
        if gatepass(rB) && (isempty(zt0) || liB{1}.drop<=zt0),  r = rB;  li = liB;  ok = true;  end
        att = sprintf('slack of the solve (%s), then the support of D one degree up',rA);
    catch ME                        % e.g. above the size gate of the harness
        att = sprintf('slack of the solve (support of D one degree up: %s)',strtok(ME.message,newline));
    end
end
pre = sprintf('%s%s; slack entries %s (solve: %d)',pre,att,mat2str(cellfun(@(s) s.ndec,li)),aux.li.ndec);
end

function m = maxabs_opvar(D)
% M = MAXABS_OPVAR(D) the largest |coefficient| of the 'opvar' object D.
m = 0;
for k = 1:numel(D.C)
    B = D.C{k};
    if isa(B,'sopvar')
        for t = 1:numel(B.params),  if ~isempty(B.params{t}),  m = max(m,full(max(abs(B.params{t}(:)))));  end,  end
    elseif isnumeric(B) && ~isempty(B)
        m = max(m,full(max(abs(B(:)))));
    end
end
end

function [E,drop] = zero_unreached(E)
% [E,DROP] = ZERO_UNREACHED(E) sets to zero the constant coefficients of the
% 'dopvar' object E that carry no decision variable, the rows "0 = residue"
% of E == 0: entries of params.A whose column of params.B is zero, and a
% fixed 'sopvar' block entirely. DROP the largest |coefficient| set to zero.
% Zeroing keeps the canonical multiplier form (a zero is canonical).
drop = 0;
for k = 1:numel(E.C)
    B = E.C{k};     [i,j] = ind2sub(size(E.C),k);
    if isa(B,'sdopvar')
        hit = false;
        for t = 1:numel(B.params.A)
            a = B.params.A{t};
            if isempty(a) || ~any(a),  continue,  end
            if isempty(B.params.B{t}),  u = a~=0;
            else,                       u = reshape(~any(B.params.B{t},1),size(a)) & a~=0;
            end
            if any(u),  drop = max(drop,full(max(abs(a(u)))));  a(u) = 0;  B.params.A{t} = a;  hit = true;  end
        end
        if hit,  E.C{i,j} = [];  E.C{i,j} = B;  end     % extract, edit, reinsert (in place)
    elseif isa(B,'sopvar')
        for t = 1:numel(B.params)
            if ~isempty(B.params{t}) && any(B.params{t}(:)),  drop = max(drop,full(max(abs(B.params{t}(:)))));  end
        end
        E.C{i,j} = [];                  % no slack variable reaches the block
    end
end
end

function [ok,r,pre] = stage2h2(PIE,S,prog,info)
% PIETOOLS_H2_norm, form 'c': R, Q, W and gam fixed (getsol).
aux = info.aux;     gv = aux.obj;
Rn = getsol_lpivar(prog,aux.Rm);    Qn = getsol_lpivar(prog,aux.Qm);    Wn = getsol_lpivar(prog,aux.Wm);
dtab = prog.decvartable;    if isstring(dtab),  dtab = cellstr(dtab);  end
gn = prog.solinfo.RRx(find(strcmp(dtab,gv.dvarname{1}),1));
PIEc = exec_pie(PIE);
st = exec_sop_settings(S);  tl = exec_tools(PIEc,st);  st.sop = tl.sop;
tr = trace_rn(Wn);
if isa(tr,'dpvar'),  trn = full(tr.C(1));  else,  trn = double(tr);  end
gc = max(gn,trn);
Iw = tl.ident(PIEc.B1,2);
Dneg = [-(gc*Iw), PIEc.B1'; PIEc.B1, Qn'*PIEc.A' + PIEc.A*Qn];
Dpos = [Wn, PIEc.C1*Qn; Qn'*PIEc.C1', Rn];
% No magnitude cleaning: slack_like zeroes the residue rows of -Dneg (see
% stage2). The executive returns no slack operator, so the residual of the
% solve that bounds the residue is unknown: a nonzero residue fails
opp = struct('like',Rn,'dw',st.sop.dw_Q+info.dw,'margin',st.sop.margin,'cap',tl.capfor(Dpos));
[r,li] = slackonly(PIEc,S,{-Dneg, Dpos},{struct('like_li',aux.li), opp});
dz = li{1}.drop;
ok = gatepass(r) && dz==0;   att = sprintf('slack of the solve, residue %.2g',dz);
if ~ok
    [rB,liB] = slackonly(PIEc,S,{-Dneg, Dpos},{struct('dw',info.dw+1,'margin',0), opp});
    if gatepass(rB),  r = rB;  li = liB;  ok = true;  end
    att = 'slack of the support, one degree up';
end
r.gam = gc;
pre = sprintf('%s; gam %.8g; trace W %.8g; slack entries %s',att,gn,trn,mat2str(cellfun(@(s) s.ndec,li)));
end

function [r,li] = slackonly(PIEc,S,tg,op)
% A program of slacks only: tg{k} >= 0 with the options op{k} (lpi_ineq),
% or, when op{k} has the field like_li, the slack recorded in like_li rebuilt
% exactly (slack_like); solved.
st = exec_sop_settings(S);  tl = exec_tools(PIEc,st);  st.sop = tl.sop;
prog = lpiprogram(PIEc.vars(:,1),PIEc.vars(:,2),PIEc.dom,st);
li = cell(size(tg));
for k = 1:numel(tg)
    if isfield(op{k},'like_li'),  [prog,li{k}] = slack_like(prog,tg{k},op{k});
    else,                         [prog,~,li{k}] = lpi_ineq(prog,tg{k},op{k});
    end
end
prog = lpisolve(prog,st.sos_opts);
cl = lpi_classify(prog,[]);
r = solrec(prog);   r.exst = cl.status;     r.eta = cl.eta;
end

function [prog,li] = slack_like(prog,P,o)
% [PROG,LI] = SLACK_LIKE(PROG,P,O) the slack of the solve on the target P:
% the positive operator of lpi_ineq with the degree specification, the
% terms and the kept basis operators (li.include, its pruning) that O.like_li
% records, on the spaces of O.ref (the operator of the solve), and the same
% margin; P - Pop == 0 imposed, as lpi_ineq does, with its residue rows
% zeroed (zero_unreached; LI.drop the largest value zeroed, LI.ndec the
% number of slack variables of this rebuild).
li = o.like_li;   ref = P;   if isfield(o,'ref') && ~isempty(o.ref),  ref = o.ref;  end
[sp,dm] = opvar_space_list(ref,'out');      % the spaces of the operator of the solve
dom = struct('vars',{ref.vars},'dom',ref.dom);
popts = struct('psatz',li.terms.codes,'psatz_offset',li.terms.offsets);
if isfield(li.terms,'rows') && ~isempty(li.terms.rows)
    popts.psatz_rows = li.terms.rows;   popts.psatz_rows_offset = li.terms.rows_offsets;
end
if ~isempty(li.include) && any(~cellfun(@isempty,li.include)),  popts.include = li.include;  end
[prog,Pop] = poslpivar_direct(prog,dm,sp,dom,li.deg,popts);
if li.margin>0,  P = P + li.margin*li.eD;  end
[E,li.drop] = zero_unreached(P - Pop);
prog = lpi_eq(prog,E,'symmetric');
li.ndec = 0;                        % every 'sdopvar' block holds the one list of the container
for k = 1:numel(Pop.C),  if isa(Pop.C{k},'sdopvar'),  li.ndec = numel(Pop.C{k}.Zd);  break,  end,  end
end

function ok = gatepass(r)
% The harness gate: solver exit 'certified', executive optimal, eta <= 1e-7.
ok = strcmp(fs(r,'status'),'certified') && strcmp(fs(r,'exst'),'optimal') && fn(r,'eta') <= 1e-7;
end

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

function s = fs(r,k),  s = '';   if isstruct(r) && isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = fn(r,k),  x = NaN;  if isstruct(r) && isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end
