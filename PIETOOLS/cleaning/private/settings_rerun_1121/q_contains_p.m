function q_contains_p(key,SP,SQ,tag,variants,room)
% Q_CONTAINS_P(KEY,SP,SQ,TAG) tests where the 'Q' LPI of
% PIETOOLS_stability_sop excludes a 'P' certificate of the catalog row KEY.
% With Q = P*T, sym(T'*Q) = T'*P*T >= eppos2*T'*T whenever P >= eppos2*I,
% and A'*Q + Q'*A = T'*P*A + A'*P*T, so every 'P' certificate is a 'Q'
% certificate when the 'Q' program can represent it. The storage P of the
% 'P' solve at the settings SP is fixed (getsol), Q = P*T, and against the
% 'Q' program at the settings SQ the script checks:
%   (a) Q = P*T in the span of the decision Q (lpivar_cdopvar with the
%       degrees deg_Q of SQ): the largest coefficient of P*T that no
%       variable of Q reaches;
%   (b) T'*P*T - eppos2*T'*T >= 0 with the slack of the margin of the 'Q'
%       solve (info.aux.lip, Npm), rebuilt exactly;
%   (c) -(T'*P*A + A'*P*T) >= 0 with the slack of -D of the 'Q' solve
%       (info.aux.li), and, as control, with that of the 'P' solve.
% (b) and (c) are programs of slack variables only, gated as post_verify
% gates. Settings rerun of 11/21/2026, maintainer 11/22: "'Q' certified at
% no setting tried ... this shouldn't happen, since P is a special case of
% Q". Results: q_contains_p_<TAG>.csv in this folder.
%
% INPUT
% - key:  char, catalog key of a stability row (executive PIE2PDEstability);
% - SP:   'struct', settings of the 'P' solve;
% - SQ:   'struct', settings of the 'Q' solve;
% - tag:  char, the name of the results file;
% - variants: (optional) logical, also the checks with room in the target
%         (weak feasibility: with P fixed, a Gram matrix of the slack forced
%         to the boundary of the cone leaves no interior point): (b) at
%         eppos2/10 and at 0; (c) -D + 1e-6 max|D| e_D, e_D the lift
%         identity of the slack of the 'Q' solve;
% - room: (optional) row of c > 0: (b) and (c) again with c*s*T'*T added to
%         the target, s the ratio of its largest coefficient to that of
%         T'*T; the least c that certifies measures how far the target is
%         from a Gram of the 'Q' slack that is positive semidefinite.
% MMP, 11/22/2026.
if nargin<5,  variants = false;  end
if nargin<6,  room = [];  end
here = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/testing_codes/private/baselines/2026-11-18_studies/private';
out = fullfile(fileparts(mfilename('fullpath')),['q_contains_p_',tag,'.csv']);
cd(here);  lb_setup();  cd(here);  addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
c = Cat(find(strcmp({Cat.key},key),1));
PIE = lb_build_system(c,[],true);
fprintf('## Q_CONTAINS_P %s %s  %s\n',tag,key,char(datetime('now')));
% The 'P' solve and the 'Q' solve, kept for their operators
clearmaps();
setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
setappdata(0,'LB_KEEP_LAST',true);
rP = lb_exec(['sop:',c.exec],PIE,SP);   LP = getappdata(0,'LB_LAST');
fprintf('## P solve: %s exst %s fr %.3g eta %.2g ndec %d\n',rP.status,rP.exst,rP.fr,rP.eta,rP.ndec);
clearmaps();
setappdata(0,'LB_KEEP_LAST',true);
rQ = lb_exec(['sop:',c.exec],PIE,SQ);   LQ = getappdata(0,'LB_LAST');
setappdata(0,'LB_KEEP_LAST',false);
fprintf('## Q solve: %s exst %s fr %.3g eta %.2g ndec %d\n',rQ.status,rQ.exst,rQ.fr,rQ.eta,rQ.ndec);
PIEc = exec_pie_sop(PIE);   T = PIEc.T;     A = PIEc.A;
auxP = LP.info.aux;     auxQ = LQ.info.aux;
Pn = getsol_lpivar_sop(LP.prog,auxP.Pm);          % P of the 'P' solve (eppos I included)
Qn = Pn*T;                                          % Q = P*T
PQ = T'*Qn;                                         % sym(T'*Q) = T'*P*T
ep = LQ.info.eppos;     e2 = min(ep(min(2,end):end));
DP = A'*Qn + Qn'*A;                                 % = T'*P*A + A'*P*T
if LQ.info.epneg~=0,  DP = DP + 2*LQ.info.epneg*PQ;  end
% (a) P*T in the span of the decision Q of the 'Q' program
stQ = exec_sop_settings(LQ.S);  tl = exec_tools_sop(PIEc,stQ);  stQ.sop = tl.sop;
qdeg = stQ.sop.deg_Q;   if isempty(qdeg),  qdeg = stQ.lpi.var.deg;  end
if ~iscell(qdeg),  qdeg = {qdeg};  end
qdeg = qdeg{min(LQ.info.nv,end)};
dwq = 0;   if ~strcmp(stQ.sop.raise_Q,'slack'),  dwq = LQ.info.dw;  end
q = qdeg;  if isnumeric(q),  q(1) = q(1) + dwq;  end
prog0 = lpiprogram(PIEc.vars(:,1),PIEc.vars(:,2),PIEc.dom,stQ);
[~,Qm] = lpivar_cdopvar(prog0,T,q);
[~,dropa] = zero_unreached(Qn - Qm);   qmax = maxabs_copvar(Qn);
fprintf('## (a) deg_Q %s: largest coefficient of P*T outside the span of Q %.3g (max|P*T| %.3g)\n',mat2str(q),dropa,qmax);
% (b) the margin with the slack of the margin of the 'Q' solve
rows = {};
if isfield(auxQ,'lip') && ~isempty(auxQ.lip)
    [rb,lib] = slackonly(PIEc,LQ.S,{PQ - e2*(T'*T)},{struct('like_li',auxQ.lip,'ref',auxQ.Pm)});
    rows(end+1,:) = {'b margin, Q slack',rb,lib{1}};
end
% (c) -D with the slack of -D of the 'Q' solve; the 'P' slack as control
[rc,lic] = slackonly(PIEc,LQ.S,{-DP},{struct('like_li',auxQ.li,'ref',auxQ.Dm)});
rows(end+1,:) = {'c -D, Q slack',rc,lic{1}};
[rp,lipp] = slackonly(PIEc,LP.S,{-DP},{struct('like_li',auxP.li,'ref',auxP.Dm)});
rows(end+1,:) = {'c -D, P slack (control)',rp,lipp{1}};
if variants                         % room in the target: weak feasibility
    if isfield(auxQ,'lip') && ~isempty(auxQ.lip)
        [r1,l1] = slackonly(PIEc,LQ.S,{PQ - 0.1*e2*(T'*T)},{struct('like_li',auxQ.lip,'ref',auxQ.Pm)});
        rows(end+1,:) = {'b margin eppos2/10, Q slack',r1,l1{1}};
        [r2,l2] = slackonly(PIEc,LQ.S,{PQ},{struct('like_li',auxQ.lip,'ref',auxQ.Pm)});
        rows(end+1,:) = {'b margin 0, Q slack',r2,l2{1}};
        % controls on the same slack: a congruence by T alone, and the
        % identity (strictly positive in every multiplier direction)
        [r4,l4] = slackonly(PIEc,LQ.S,{T'*T},{struct('like_li',auxQ.lip,'ref',auxQ.Pm)});
        rows(end+1,:) = {'b control T''*T, Q slack',r4,l4{1}};
        [r5,l5] = slackonly(PIEc,LQ.S,{eye_copvar_sop(auxQ.Pm)},{struct('like_li',auxQ.lip,'ref',auxQ.Pm)});
        rows(end+1,:) = {'b control I, Q slack',r5,l5{1}};
    end
    % the same targets with slacks sized from their own (numeric) supports,
    % as lpi_ineq_sop sizes them (pruned to the nonzero cells), raise 0:
    % the Newton polytope of the 'P' point, not of the generic decision Q
    [r6,l6] = slackonly(PIEc,LQ.S,{PQ - e2*(T'*T)},{struct('dw',0,'margin',0)});
    rows(end+1,:) = {'b margin, slack from its own support',r6,addrop(l6{1})};
    [r7,l7] = slackonly(PIEc,LQ.S,{-DP},{struct('dw',0,'margin',0)});
    rows(end+1,:) = {'c -D, slack from its own support',r7,addrop(l7{1})};
    % -D with a slack sized like the generalized storage sym(T'*Q) of the
    % 'Q' program (opts.like, as the 'P' executive sizes its slack like P),
    % weight raises 0, 1, 2
    for dwl = 0:2
        [r8,l8] = slackonly(PIEc,LQ.S,{-DP},{struct('like',auxQ.Pm,'dw',dwl,'margin',0)});
        rows(end+1,:) = {sprintf('c -D, slack like sym(T''*Q), dw %d',dwl),r8,addrop(l8{1})}; %#ok<AGROW>
    end
    % the margin with a slack sized like sym(T'*Q), weight raises 0, 1, 2
    for dwl = 0:2
        [r9,l9] = slackonly(PIEc,LQ.S,{PQ - e2*(T'*T)},{struct('like',auxQ.Pm,'dw',dwl,'margin',0)});
        rows(end+1,:) = {sprintf('b margin, slack like sym(T''*Q), dw %d',dwl),r9,addrop(l9{1})}; %#ok<AGROW>
    end
    if isfield(auxQ.li,'eD') && ~isempty(auxQ.li.eD)
        dmax = maxabs_copvar(DP);
        [r3,l3] = slackonly(PIEc,LQ.S,{-DP + (1e-6*dmax)*auxQ.li.eD},{struct('like_li',auxQ.li,'ref',auxQ.Dm)});
        rows(end+1,:) = {'c -D + 1e-6 max|D| eD, Q slack',r3,l3{1}};
    end
end
if ~isempty(room)                   % room c*s*T'*T in the target
    TT = T'*T;  tt = maxabs_copvar(TT);
    sb = maxabs_copvar(PQ)/tt;      sc = maxabs_copvar(DP)/tt;
    for cr = room
        if isfield(auxQ,'lip') && ~isempty(auxQ.lip)
            [rr,ll] = slackonly(PIEc,LQ.S,{PQ - e2*TT + (cr*sb)*TT},{struct('like_li',auxQ.lip,'ref',auxQ.Pm)});
            rows(end+1,:) = {sprintf('b margin + %g s T''*T, Q slack',cr),rr,ll{1}}; %#ok<AGROW>
        end
        [rr,ll] = slackonly(PIEc,LQ.S,{-DP + (cr*sc)*TT},{struct('like_li',auxQ.li,'ref',auxQ.Dm)});
        rows(end+1,:) = {sprintf('c -D + %g s T''*T, Q slack',cr),rr,ll{1}}; %#ok<AGROW>
    end
end
fid = fopen(out,'w');
fprintf(fid,'check,status,fr,exst,eta,ndec,residue_zeroed,slack_ndec\n');
fprintf(fid,'a span of Q (deg_Q %s),,,,,,%.3g,\n',strrep(mat2str(q),' ',' '),dropa);
for i = 1:size(rows,1)
    r = rows{i,2};  li = rows{i,3};
    fprintf(fid,'%s,%s,%.4g,%s,%.3g,%g,%.3g,%d\n',rows{i,1},r.status,r.fr,r.exst,r.eta,r.ndec,li.drop,li.ndec);
    fprintf('## (%s): %s fr %.3g exst %s eta %.2g; residue zeroed %.2g; slack %d\n',rows{i,1},r.status,r.fr,r.exst,r.eta,li.drop,li.ndec);
end
fclose(fid);
fprintf('## Q_CONTAINS_P %s DONE %s\n',tag,char(datetime('now')));
end

function li = addrop(li)
% The info of lpi_ineq_sop with the fields the report prints (no residue
% is zeroed on that path; ndec as lpi_ineq_sop counts it).
li.drop = 0;
end

function clearmaps()
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_LAST','LB_KEEP_LAST'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
end

% ------------------------------------------------------------------------
% The slack programs of post_verify (2026-11-18_studies harness, 11/21)
function [r,li] = slackonly(PIEc,S,tg,op)
st = exec_sop_settings(S);  tl = exec_tools_sop(PIEc,st);  st.sop = tl.sop;
prog = lpiprogram(PIEc.vars(:,1),PIEc.vars(:,2),PIEc.dom,st);
li = cell(size(tg));
for k = 1:numel(tg)
    if isfield(op{k},'like_li'),  [prog,li{k}] = slack_like(prog,tg{k},op{k});
    else,                         [prog,~,li{k}] = lpi_ineq_sop(prog,tg{k},op{k});
    end
end
prog = lpisolve(prog,st.sos_opts);
cl = lpi_classify_sop(prog,[]);
r = solrec(prog);   r.exst = cl.status;     r.eta = cl.eta;
end

function [prog,li] = slack_like(prog,P,o)
li = o.like_li;   ref = P;   if isfield(o,'ref') && ~isempty(o.ref),  ref = o.ref;  end
[sp,dm] = copvar_space_list(ref,'out');
dom = struct('vars',{ref.vars},'dom',ref.dom);
popts = struct('psatz',li.terms.codes,'psatz_offset',li.terms.offsets);
if isfield(li.terms,'rows') && ~isempty(li.terms.rows)
    popts.psatz_rows = li.terms.rows;   popts.psatz_rows_offset = li.terms.rows_offsets;
end
if ~isempty(li.include) && any(~cellfun(@isempty,li.include)),  popts.include = li.include;  end
[prog,Pop] = poscopvar_direct(prog,dm,sp,dom,li.deg,popts);
if li.margin>0,  P = P + li.margin*li.eD;  end
[E,li.drop] = zero_unreached(P - Pop);
prog = lpi_eq_sop(prog,E,'symmetric');
li.ndec = 0;
for k = 1:numel(Pop.C),  if isa(Pop.C{k},'sdopvar'),  li.ndec = numel(Pop.C{k}.Zd);  break,  end,  end
end

function m = maxabs_copvar(D)
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
% The constants of E == 0 that carry no decision variable set to zero; DROP
% the largest of them.
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
        if hit,  E.C{i,j} = [];  E.C{i,j} = B;  end
    elseif isa(B,'sopvar')
        for t = 1:numel(B.params)
            if ~isempty(B.params{t}) && any(B.params{t}(:)),  drop = max(drop,full(max(abs(B.params{t}(:)))));  end
        end
        E.C{i,j} = [];
    end
end
end

function r = solrec(prog)
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
