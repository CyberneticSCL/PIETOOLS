function R = test_getsol_sop(parts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_GETSOL_SOP(PARTS) tests getsol_lpivar_sop, lpigetsol_sop and
% subs_dvar_sop against the SEMANTICS of the classes, never against an
% inverse routine (CLAUDE.md sec. 4). Errors on the first failure; returns
% the measured numbers. PARTS is a cellstr subset of {'i','ii','iii','iv'}
% (default all).
%
% (i)   Random 'sdopvar'/'cdopvar' at a random decision vector: 1-D, 2-D and
%       3 variables, R^n blocks, structurally empty blocks, per-block name
%       lists merged by the constructor, a 'sopvar' block among 'sdopvar'
%       ones. The name list is shuffled and padded with names absent from
%       the operator, and RRx carries entries past decvartable (the slack
%       a solve appends). Each extracted block is applied to random
%       polynomial test functions by quadrature straight from the 'sopvar'
%       class definition (heatNd_apply) and compared with an INDEPENDENT
%       evaluation: the block's gamma cells converted by 'sdvar2dpvar' (DJ)
%       to 'dpvar' kernels, the values substituted by legacy 'sosgetsol',
%       and the polynomial kernels integrated by a separate quadrature.
%       A control with the values shuffled among the names must fail the
%       same comparison. Also: partial substitution (names left out stay
%       decision variables, in order, with their B rows), class invariants
%       (verify, no canonical-form rewrite), and agreement of the three
%       entry points.
% (ii)  A real solved program: cx_PIE2PDEstability on cx_plant('rd',0.5),
%       MOSEK. Certified first (MOSEK status, rel_b on solinfo.RRx rows,
%       PSD Gram blocks). Every decision operator (P, Q, D, N) is
%       extracted; both lpi_eq relations, T'*Q - P = 0 and D + N = 0, are
%       evaluated on the extracted operators by their action on random
%       test functions, relative to the size of the two sides; P and N
%       must be PSD on test functions. Controls: solinfo.x in place of
%       RRx, and a perturbed RRx, must both break the relations.
% (iii) Dispatch: every legacy class through lpigetsol_sop against legacy
%       lpigetsol on a program with a known decision vector (bit-identical
%       results), the same for subs_dvar_sop, the new classes, a cell
%       mixing both families, and errors for unknown classes and for an
%       unsolved program; a name listed twice in decvartable reads its
%       lowest row, as getequation places its At row and sosgetsol reads.
% (iv)  A solved program with a free scalar and an inequality, where
%       solinfo.x and RRx differ in length: the 1-D primal H-infinity LPI
%       on cx_plant('io1') with gamma a decision variable (gam >= 0,
%       minimized), MOSEK. The extracted KYP operator must equal the one
%       rebuilt from the legacy gamma (lpigetsol on the dpvar) and the
%       extracted Q, be negative on test functions, and T'*Q = R must hold
%       on the extracted operators. Uses the dpvar-operator overloads
%       (-gam*Iw, [ ]) and eye_copvar_sop of Tier 1b/1c.
%
% Requires the PIETOOLS path (pietools_path_update) and, for (ii) and
% (iv), MOSEK and the cx_exec harness
% (sopvar/Testfolder/sdopvar/claude_tests/cx_exec).
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(parts),  parts = {'i','ii','iii','iv'};  end
if ischar(parts),   parts = {parts};    end
R = struct();
if any(strcmp(parts,'i')),      rng(20260929);  R.i   = part_i();     end
if any(strcmp(parts,'ii')),     rng(20260930);  R.ii  = part_ii();    end
if any(strcmp(parts,'iii')),    rng(20261001);  R.iii = part_iii();   end
if any(strcmp(parts,'iv')),     rng(20261002);  R.iv  = part_iv();    end
fprintf('test_getsol_sop: parts {%s} passed.\n',strjoin(parts,','));
end


% ========================================================================
% (i) random operators against an independent evaluation
% ========================================================================
function out = part_i()
tol = 1e-10;    nq = 6;     % Gauss nodes: exact to degree 11, integrands here <= 4
V = {'s1','s2','s3'};   D = [0 1; -1 2; 0.5 1.5];   % distinct domains per variable
sp = @(o,i) struct('out',{o},'in',{i},'vars',{V});
% {label, spaces, dims, deg, ndec, occ, sopvar block (i,j) or []}
cases = {
  {'1-D s1->s1',            sp({{'s1'}},{{'s1'}}),                     struct('out',2,'in',3),         2, 9,       [],  []}
  {'1-D s2->s1 (S1,S2)',    sp({{'s1'}},{{'s2'}}),                     struct('out',2,'in',1),         2, 5,       [],  []}
  {'2-D s1s2->s1s2',        sp({{'s1','s2'}},{{'s1','s2'}}),           struct('out',1,'in',2),         1, 6,       [],  []}
  {'s2s3->s1s2',            sp({{'s1','s2'}},{{'s2','s3'}}),           struct('out',2,'in',1),         1, 6,       [],  []}
  {'3 vars s1s2s3',         sp({{'s1','s2','s3'}},{{'s1','s2','s3'}}), struct('out',1,'in',1),         1, 5,       [],  []}
  {'R^n x L2 2x2',          sp({{},{'s1'}},{{},{'s1'}}),               struct('out',[2;1],'in',[3;1]), 2, 7,       [],  []}
  {'3x3 empty blocks',      sp({{},{'s1'},{'s1','s2'}},{{'s1'},{},{'s2','s3'}}), ...
                            struct('out',[2;1;1],'in',[1;2;1]),        1, 8, logical([1 0 1;1 1 0;0 1 1]), []}
  {'2x2 per-block lists',   sp({{'s1'},{}},{{'s1'},{'s2'}}),           struct('out',[1;2],'in',[2;1]), 2, [3 4;5 2],[],  []}
  {'2x2 sopvar block',      sp({{'s1'},{}},{{'s1'},{'s2'}}),           struct('out',[1;2],'in',[2;1]), 2, 6,       [],  [1 2]}
  };
out = struct('label',{},'q',{},'nblk',{},'err',{},'err_ctrl',{},'err_part',{});
for c = 1:numel(cases)
    cs = cases{c};
    [lbl,spaces,dims,deg,ndec,occ,sblk] = cs{:};
    P = rand_cdopvar(spaces,dims,D,deg,ndec,0.5,occ);
    if ~isempty(sblk)                   % a fixed block among decision blocks
        F = rand_copvar(spaces,dims,D,deg,0.5,occ);
        C = P.C;    C{sblk(1),sblk(2)} = F.C{sblk(1),sblk(2)};
        P = cdopvar(C);
        assert(isa(P.C{sblk(1),sblk(2)},'sopvar'),'case %s: sopvar block not kept',lbl)
    end
    assert(verify(P).true,'case %s: generator produced an invalid cdopvar',lbl)
    Zd = cellstr(P.Zd(:));      q = numel(Zd);
    d = randn(q,1);
    % Shuffled list, padded with names the operator does not have.
    extra = strcat('zz_absent_',cellstr(string((1:4)')));
    names = [Zd; extra];        vals = [d; randn(4,1)];
    pm = randperm(numel(names));    names = names(pm);  vals = vals(pm);
    prog = struct('decvartable',{names},'solinfo',struct('RRx',[vals; randn(3,1)],'info',struct('pinf',0)));

    lastwarn('');
    S1 = subs_dvar_sop(P,names,vals);
    S2 = getsol_lpivar_sop(prog,P);
    S3 = lpigetsol_sop(prog,P);
    [wm,wid] = lastwarn;
    assert(~contains(wid,'noncanonical'),'case %s: result was rewritten (%s)',lbl,wm)
    assert(isa(S1,'copvar') && isa(S2,'copvar') && isa(S3,'copvar'),'case %s: not copvar',lbl)
    assert(verify(S1).true,'case %s: result fails verify',lbl)
    assert(same_ops(S1,S2,0) && same_ops(S1,S3,0),'case %s: entry points disagree',lbl)
    assert(isequal(S1.vars,P.vars) && isequal(S1.dom,P.dom) && isequal(S1.space_out,P.space_out) ...
        && isequal(S1.space_in,P.space_in) && isequal(S1.dim_out,P.dim_out) && isequal(S1.dim_in,P.dim_in), ...
        'case %s: container metadata changed',lbl)

    % Control: the same values shuffled among the operator names.
    dperm = d(derangement(q));
    Sbad = subs_dvar_sop(P,Zd,dperm);

    [M,N] = size(P.C);  err = 0;    errc = inf;     nblk = 0;
    for i = 1:M
        for j = 1:N
            b = P.C{i,j};   bs = S1.C{i,j};
            if isempty(b)
                assert(isempty(bs),'case %s: zero block (%d,%d) not kept',lbl,i,j);    continue
            end
            if isa(b,'sopvar')
                assert(isequal(bs,b),'case %s: sopvar block (%d,%d) changed',lbl,i,j); continue
            end
            assert(isa(bs,'sopvar'),'case %s: block (%d,%d) is %s',lbl,i,j,class(bs))
            % the standalone sdopvar path gives the same block
            Sb = subs_dvar_sop(b,names,vals);
            assert(isa(Sb,'sopvar') && same_ops(Sb,bs,0),'case %s: sdopvar path differs at (%d,%d)',lbl,i,j)
            [f,S] = test_data(b,D,V);
            y_get = heatNd_apply(bs,f,S,nq);
            y_ind = indep_apply(b,names,vals,f,S,nq);
            y_bad = heatNd_apply(Sbad.C{i,j},f,S,nq);
            nr = norm(y_ind,'fro');
            assert(nr>1e-6,'case %s: vacuous comparison at (%d,%d)',lbl,i,j)
            err = max(err,norm(y_get-y_ind,'fro')/nr);
            errc = min(errc,norm(y_bad-y_ind,'fro')/nr);
            nblk = nblk+1;
        end
    end
    assert(err<tol,'case %s: extracted operator differs from the independent evaluation by %.3g',lbl,err)
    assert(errc>1e-3,'case %s: control not detected (%.3g) - the comparison cannot see a wrong d',lbl,errc)

    % Partial substitution: a random half of the names given.
    giv = randperm(q,floor(q/2));   rest = setdiff(1:q,giv);
    Sp = subs_dvar_sop(P,Zd(giv),d(giv));
    assert(isa(Sp,'cdopvar') && verify(Sp).true,'case %s: partial result invalid',lbl)
    assert(isequal(cellstr(Sp.Zd(:)),Zd(rest)),'case %s: remaining names wrong or reordered',lbl)
    ep = 0;
    for k = 1:numel(P.C)
        if ~isa(P.C{k},'sdopvar'),  continue,   end
        Bo = P.C{k}.params.B;   Bp = Sp.C{k}.params.B;
        for g = 1:numel(Bo)
            ep = max(ep,norm(Bp{g}-Bo{g}(rest,:),1));   % kept rows are the original rows
        end
    end
    Sf = subs_dvar_sop(Sp,Zd,d);                        % then the rest
    assert(isa(Sf,'copvar') && same_ops(Sf,S1,1e-13),'case %s: two-step substitution differs',lbl)
    assert(ep==0,'case %s: kept B rows changed',lbl)

    out(end+1) = struct('label',lbl,'q',q,'nblk',nblk,'err',err,'err_ctrl',errc,'err_part',ep); %#ok<AGROW>
    fprintf('(i) %-22s q=%3d blocks=%d  rel err %.2e  control %.2e\n',lbl,q,nblk,err,errc);
end
end


% ------------------------------------------------------------------------
function [f,S] = test_data(b,D,V)
% Random polynomial test function on the block's input space (handle of
% points whose columns are in b.vars.in order) and three evaluation points
% of its output space (columns in b.vars.out order), by variable NAME.
vin = b.vars.in(:)';    vout = b.vars.out(:)';
p = b.dims(2);
nt = 4;                                         % terms, degree <= 2 per variable
F.deg = randi([0 2],nt,numel(vin));     F.deg(1,:) = 0;
F.coef = randn(nt,p);
f = @(X) prod_terms(X,F);
S = zeros(3,numel(vout));
for k = 1:numel(vout)
    dk = D(strcmp(V,vout{k}),:);
    S(:,k) = dk(1) + (dk(2)-dk(1))*rand(3,1);
end
if isempty(vout),   S = zeros(1,0);     end
end

function y = prod_terms(X,F)
n = size(X,1);  Z = ones(n,size(F.deg,1));
for t = 1:size(F.deg,1)
    for k = 1:size(F.deg,2),    Z(:,t) = Z(:,t).*X(:,k).^F.deg(t,k);    end
end
y = Z*F.coef;
end


% ------------------------------------------------------------------------
function y = indep_apply(b,names,vals,f,S,nq)
% (P(d) x)(s) for an 'sdopvar' block b, with its kernels obtained without
% the code under test: sdvar2dpvar gives each gamma cell as a 'dpvar' in
% (vars.out, vars.in_dum); legacy sosgetsol substitutes the values by
% name; the polynomial kernels are integrated here. Gamma cell g over the
% SORTED shared variables, direction 1 fastest; 1 = multiplier (s' = s),
% 2 = int_a^s, 3 = int_s^b (sopvar.m header, sopvar document Sec. 4).
vin = b.vars.in(:)';    vout = b.vars.out(:)';
vdum = strcat(vin,'_dum');
prm = b.params;     Zd = cellstr(b.Zd(:));
sos = struct('decvartable',{names(:)},'solinfo',struct('RRx',vals(:)));
S3 = sort(intersect(vin,vout));     n3 = numel(S3);
[~,i3in] = ismember(S3,vin);        [~,i3out] = ismember(S3,vout);
S1 = setdiff(vin,vout);             [~,i1in] = ismember(S1,vin);
q = b.dims(1);  p = b.dims(2);
[xg,wg] = gauss01(nq);
ns = size(S,1);     y = zeros(ns,q);
for g = 1:numel(prm.A)
    Dg = sdvar2dpvar(struct('A',prm.A{g},'B',prm.B{g}),b.dims, ...
                     struct('out',{vout},'in',{vdum}),b.ZL,b.ZR,Zd);
    K = sosgetsol(sos,Dg);                      % q x p, polynomial or double
    gam = ones(1,n3);
    if n3>0,    c = cell(1,n3);     [c{:}] = ind2sub([3*ones(1,n3),1],g);   gam = [c{1:n3}];    end
    for is = 1:ns
        s = S(is,:);
        lo = [];    hi = [];    col = [];
        for k = 1:n3
            if gam(k)==2,       lo(end+1) = b.dom.in(i3in(k),1); hi(end+1) = s(i3out(k));        col(end+1) = i3in(k); %#ok<AGROW>
            elseif gam(k)==3,   lo(end+1) = s(i3out(k));         hi(end+1) = b.dom.in(i3in(k),2); col(end+1) = i3in(k); %#ok<AGROW>
            end
        end
        for k = 1:numel(S1)
            lo(end+1) = b.dom.in(i1in(k),1);   hi(end+1) = b.dom.in(i1in(k),2);   col(end+1) = i1in(k); %#ok<AGROW>
        end
        nd = numel(col);    nx = nq^nd;
        X = zeros(nx,numel(vin));   w = ones(nx,1);
        for k = 1:n3
            if gam(k)==1,   X(:,i3in(k)) = s(i3out(k));   end
        end
        sub = cell(1,nd);
        if nd>0,    [sub{:}] = ind2sub([nq*ones(1,nd),1],(1:nx)');  end
        for t = 1:nd
            X(:,col(t)) = lo(t) + (hi(t)-lo(t))*xg(sub{t});
            w = w.*((hi(t)-lo(t))*wg(sub{t}));
        end
        Kv = poly_at(K,[vout,vdum],[repmat(s,nx,1),X]);    % nx x q*p, column-major
        fx = f(X);
        for ii = 1:q
            for jj = 1:p
                y(is,ii) = y(is,ii) + sum(w.*Kv(:,ii+q*(jj-1)).*fx(:,jj));
            end
        end
    end
end
end

function v = poly_at(K,names,Y)
% Values of a matrix polynomial (or double) at the rows of Y, whose columns
% are the variables NAMES; one column per entry, column-major.
n = size(Y,1);
if isa(K,'double')
    v = repmat(reshape(full(K),1,[]),n,1);  return
end
dm = full(K.degmat);    cf = full(K.coefficient);   vn = K.varname;
[tf,loc] = ismember(vn,names);
assert(all(tf),'poly_at: kernel variable not among %s',strjoin(names,','))
Z = ones(n,size(dm,1));
for t = 1:size(dm,1)
    for k = 1:numel(loc),   Z(:,t) = Z(:,t).*Y(:,loc(k)).^dm(t,k);  end
end
v = Z*cf;
end


% ========================================================================
% (ii) a real solved container stability program
% ========================================================================
function out = part_ii()
nqi = 10;   nqo = 12;   ntest = 6;
PIE = cx_plant('rd',0.5);
st = cx_settings('light','mosek');
t0 = tic;   [prog,ops] = cx_PIE2PDEstability(PIE,st);    tb = toc(t0);
r = cx_solve(prog,st.sos_opts);
sol = r.sol;
fprintf('(ii) rd 0.5: build %.1f s, solve %.1f s wall, pinf %g dinf %g numerr %g, rel_b %.2e, psd_min %.2e (rel %.2e)\n',...
    tb,r.wall,r.pinf,r.dinf,r.numerr,r.rel_b,r.psd_min,r.psd_relmin);
assert(r.st==1,'(ii) MOSEK did not certify feasibility (pinf %g, numerr %g)',r.pinf,r.numerr)
assert(r.rel_b<=1e-6,'(ii) equality residual rel_b = %.3g on solinfo.RRx rows',r.rel_b)
assert(r.psd_relmin>=-1e-8,'(ii) Gram blocks not PSD: min eig / max eig = %.3g',r.psd_relmin)
ntot = numel(sol.decvartable);
nxr = min(ntot,numel(sol.solinfo.x));
fprintf('(ii) %d decision variables, RRx %d, x %d; ||x-RRx|| = %.3g on the first %d, ||RRx|| = %.3g\n',ntot,...
    numel(sol.solinfo.RRx),numel(sol.solinfo.x),norm(sol.solinfo.x(1:nxr)-sol.solinfo.RRx(1:nxr)),nxr,...
    norm(sol.solinfo.RRx(1:ntot)));

% Extract every decision operator, through all three entry points.
nm = {'P','Q','D','N'};     X = struct();
t0 = tic;
for k = 1:numel(nm)
    Xk = getsol_lpivar_sop(sol,ops.(nm{k}));
    assert(isa(Xk,'copvar') && verify(Xk).true,'(ii) %s: extracted operator invalid',nm{k})
    assert(same_ops(Xk,lpigetsol_sop(sol,ops.(nm{k})),0),'(ii) %s: lpigetsol_sop differs',nm{k})
    assert(same_ops(Xk,subs_dvar_sop(ops.(nm{k}),sol.decvartable,sol.solinfo.RRx(1:ntot)),0),...
        '(ii) %s: subs_dvar_sop differs',nm{k})
    X.(nm{k}) = Xk;
end
tx = toc(t0);
T = ops.T;  A = ops.A;
epn = st.epneg;
% Test functions on T's input space (P, Q, N all live there).
F = test_funcs(X.P,ntest);
G = out_grid(X.P,nqo);
% Relation 1: T'*Q - P = 0 (lpi_eq_cdopvar, not symmetric).
L1 = apply_c(T'*X.Q,F,G,nqi);    R1 = apply_c(X.P,F,G,nqi);
% Relation 2: D + N = 0, D = A'*Q + Q'*A + epneg*P ('symmetric').
Dq = A'*X.Q + X.Q'*A;   if epn~=0,  Dq = Dq + epn*X.P;  end
L2 = apply_c(X.D,F,G,nqi);      R2 = apply_c(X.N,F,G,nqi);    L2q = apply_c(Dq,F,G,nqi);
rel1 = rel_res(L1,R1,G);        rel2 = rel_res(L2,neg(R2),G);
relD = rel_res(L2,L2q,G);       % getsol(D) against D rebuilt from getsol(Q)
% PSD on test functions: <x, P x> and <x, N x>, relative to ||x|| ||P x||.
pP = min_quad(F,R1,G);      pN = min_quad(F,R2,G);
fprintf('(ii) extract %.2f s | T''Q-P: rel %.2e | D+N: rel %.2e | D vs A''Q+Q''A: rel %.2e | min <x,Px> %.2e, <x,Nx> %.2e\n',...
    tx,rel1,rel2,relD,pP,pN);
assert(rel1<=1e-6 && rel2<=1e-6,'(ii) extracted operators violate the equality constraints (%.3g, %.3g)',rel1,rel2)
assert(relD<=1e-10,'(ii) getsol(D) differs from A''*getsol(Q)+getsol(Q)''*A by %.3g',relD)
assert(pP>=-1e-8 && pN>=-1e-8,'(ii) extracted P or N not PSD on test functions (%.3g, %.3g)',pP,pN)

% Controls: each must break the relations, else the check above is blind.
ctl = struct('name',{'solinfo.x','RRx+1e-3'},'rel1',{NaN,NaN},'rel2',{NaN,NaN});
xx = zeros(ntot,1);     nx = min(ntot,numel(sol.solinfo.x));    xx(1:nx) = sol.solinfo.x(1:nx);
bad = sol;  bad.solinfo.RRx = xx;
bad2 = sol; bad2.solinfo.RRx = sol.solinfo.RRx + 1e-3*norm(sol.solinfo.RRx)/sqrt(ntot)*randn(size(sol.solinfo.RRx));
bb = {bad,bad2};
for k = 1:2
    Qb = getsol_lpivar_sop(bb{k},ops.Q);    Pb = getsol_lpivar_sop(bb{k},ops.P);
    Nb = getsol_lpivar_sop(bb{k},ops.N);    Db = getsol_lpivar_sop(bb{k},ops.D);
    ctl(k).rel1 = rel_res(apply_c(T'*Qb,F,G,nqi),apply_c(Pb,F,G,nqi),G);
    ctl(k).rel2 = rel_res(apply_c(Db,F,G,nqi),neg(apply_c(Nb,F,G,nqi)),G);
    fprintf('(ii) control %-9s: T''Q-P rel %.2e, D+N rel %.2e\n',ctl(k).name,ctl(k).rel1,ctl(k).rel2);
end
assert(max(ctl(2).rel1,ctl(2).rel2)>1e3*max(rel1,rel2),'(ii) perturbed RRx not detected')
out = struct('ndec',ntot,'build',tb,'solve',r.wall,'rel_b',r.rel_b,'psd_relmin',r.psd_relmin,...
    'extract',tx,'rel1',rel1,'rel2',rel2,'relD',relD,'psdP',pP,'psdN',pN,'controls',ctl);
end


% ------------------------------------------------------------------------
function F = test_funcs(P,nt)
% nt random polynomial test functions on the column spaces of container P:
% F{t}{j} = struct(vars, deg, coef), degree <= 3 per variable.
[~,N] = size(P.C);
F = cell(1,nt);
for t = 1:nt
    F{t} = cell(1,N);
    for j = 1:N
        vj = P.vars(P.space_in(j,:));
        F{t}{j} = struct('vars',{vj},'deg',randi([0 3],5,numel(vj)),'coef',randn(5,P.dim_in(j)));
    end
end
end

function G = out_grid(P,nq)
% Gauss grid (points and weights) on each row space of P; R^n: one point.
[M,~] = size(P.C);
[xg,wg] = gauss01(nq);
G = cell(1,M);
for i = 1:M
    vi = P.vars(P.space_out(i,:));      di = P.dom(P.space_out(i,:),:);
    nd = numel(vi);     nx = nq^nd;
    Xp = zeros(nx,nd);  w = ones(nx,1);
    sub = cell(1,nd);
    if nd>0,    [sub{:}] = ind2sub([nq*ones(1,nd),1],(1:nx)');  end
    for t = 1:nd
        Xp(:,t) = di(t,1) + (di(t,2)-di(t,1))*xg(sub{t});
        w = w.*((di(t,2)-di(t,1))*wg(sub{t}));
    end
    G{i} = struct('vars',{vi},'X',Xp,'w',w);
end
end

function Y = apply_c(P,F,G,nq)
% Y{t}{i} = ((P x_t)_i at G{i}.X) = sum_j heatNd_apply(P.C{i,j}, x_t,j).
[M,N] = size(P.C);
Y = cell(1,numel(F));
for t = 1:numel(F)
    Y{t} = cell(1,M);
    for i = 1:M
        Y{t}{i} = zeros(size(G{i}.X,1),P.dim_out(i));
        for j = 1:N
            b = P.C{i,j};
            if isempty(b),  continue,   end
            [~,ci] = ismember(b.vars.in,F{t}{j}.vars);      % f columns by name
            [~,co] = ismember(b.vars.out,G{i}.vars);
            f = @(Z) eval_f(Z,ci,F{t}{j});
            Y{t}{i} = Y{t}{i} + heatNd_apply(b,f,G{i}.X(:,co),nq);
        end
    end
end
end

function y = eval_f(Z,ci,Fj)
% Test function Fj at points Z whose column k is variable Fj.vars{ci(k)}.
X = zeros(size(Z,1),numel(Fj.vars));    X(:,ci) = Z;
y = prod_terms(X,Fj);
end

function r = rel_res(L,Rt,G)
% max_t ||L_t - R_t|| / max_t max(||L_t||,||R_t||), L2 norms on the grids.
num = 0;    den = 0;
for t = 1:numel(L)
    num = max(num,l2n(cellfun(@minus,L{t},Rt{t},'uni',0),G));
    den = max([den,l2n(L{t},G),l2n(Rt{t},G)]);
end
r = num/max(den,realmin);
end

function Y = neg(Y)
for t = 1:numel(Y),     Y{t} = cellfun(@uminus,Y{t},'uni',0);  end
end

function n = l2n(Y,G)
n = 0;
for i = 1:numel(Y),     n = n + sum(G{i}.w.*sum(Y{i}.^2,2));    end
n = sqrt(n);
end

function m = min_quad(F,Y,G)
% min_t <x_t, P x_t> / (||x_t|| ||P x_t||), x_t evaluated on the grids of
% the same (square) spaces.
m = inf;
for t = 1:numel(F)
    ip = 0;     xs = cell(1,numel(G));
    for i = 1:numel(G)
        [~,ci] = ismember(G{i}.vars,F{t}{i}.vars);
        xs{i} = eval_f(G{i}.X,ci,F{t}{i});
        ip = ip + sum(G{i}.w.*sum(xs{i}.*Y{t}{i},2));
    end
    m = min(m,ip/max(l2n(xs,G)*l2n(Y{t},G),realmin));
end
end


% ========================================================================
% (iv) a solved program with a free scalar and an inequality: H-infinity
% ========================================================================
function out = part_iv()
% The 1-D primal H-infinity LPI with gamma a decision variable, built as in
% test_dpvar_hinf_chain_sop (cx_Hinf_gain with gam = dpvar, gam >= 0,
% minimized). Here solinfo.x holds one more entry than decvartable (the
% slack of the 'ineq'), which is the case where x and RRx differ.
nqi = 10;   nqo = 12;   ntest = 6;
st = cx_settings('light','mosek');
PIE = initialize(cx_plant('io1'));
Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.Bw);  Cz = opvar2copvar(PIE.Cz);  Dzw = opvar2copvar(PIE.Dzw);
[spw,dmw] = cx_space_list(Bw,'in');     [spz,dmz] = cx_space_list(Cz,'out');
Iw = eye_copvar_sop(dmw,spw,PIE.dom);   Iz = eye_copvar_sop(dmz,spz,PIE.dom);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
gam = dpvar('gam');
prog = lpidecvar(prog,gam);     prog = lpi_ineq(prog,gam);  prog = lpisetobj(prog,gam);
[prog,Rm] = cx_hinf_lf(prog,Tm,PIE,st);
[sp,dm] = cx_space_list(Tm,'out');
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,cx_hinf_qdeg(Rm));
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);
Km = [-gam*Iw, Dzw', Bw'*Qm;  Dzw, -gam*Iz, Cz;  Qm'*Bw, Cz', Am'*Qm + Qm'*Am];
prog = cx_hinf_slack(prog,Km,PIE,st);
r = cx_solve(prog,st.sos_opts);     sol = r.sol;
ntot = numel(sol.decvartable);
fprintf('(iv) io1 H-inf: pinf %g numerr %g, rel_b %.2e, psd_relmin %.2e; decvartable %d, RRx %d, x %d\n',...
    r.pinf,r.numerr,r.rel_b,r.psd_relmin,ntot,numel(sol.solinfo.RRx),numel(sol.solinfo.x));
assert(r.numerr==0 && r.pinf==0 && r.rel_b<=1e-6 && r.psd_relmin>=-1e-8,'(iv) solve not certified')
g = double(lpigetsol(sol,gam));             % legacy dpvar path
Ks = getsol_lpivar_sop(sol,Km);     Qs = getsol_lpivar_sop(sol,Qm);     Rs = getsol_lpivar_sop(sol,Rm);
assert(isa(Ks,'copvar') && isa(Qs,'copvar') && isa(Rs,'copvar'),'(iv) not fixed')
% Km at the solution, rebuilt from the legacy gamma and the extracted Q.
Kf = [(-g)*Iw, Dzw', Bw'*Qs;  Dzw, (-g)*Iz, Cz;  Qs'*Bw, Cz', Am'*Qs + Qs'*Am];
F = test_funcs(Ks,ntest);   G = out_grid(Ks,nqo);
YK = apply_c(Ks,F,G,nqi);
relK = rel_res(YK,apply_c(Kf,F,G,nqi),G);
% T'*Q - R = 0 on the extracted operators; KYP operator <= 0.
FT = test_funcs(Rs,ntest);  GT = out_grid(Rs,nqo);
rel1 = rel_res(apply_c(Tm'*Qs,FT,GT,nqi),apply_c(Rs,FT,GT,nqi),GT);
pK = -min_quad(F,neg(YK),G);                % max <x,Kx>/(||x|| ||Kx||)
fprintf('(iv) gamma %.8g | getsol(Km) vs Km(gamma, getsol(Q)): rel %.2e | T''Q-R: rel %.2e | max <x,Kx> %.2e\n',...
    g,relK,rel1,pK);
assert(relK<=1e-10,'(iv) getsol(Km) differs from Km rebuilt with the legacy gamma (%.3g)',relK)
% 1e-5, not (ii)'s 1e-6: this solve stops at rel_b ~5e-7 (measured 09/29),
% and the operator residual measured 8.4e-7 sits at that level.
assert(rel1<=1e-5,'(iv) extracted operators violate T''Q = R (%.3g)',rel1)
assert(pK<=1e-8,'(iv) extracted KYP operator not negative on test functions (%.3g)',pK)
% Control: solinfo.x in RRx's place (x is the solver's cone vector).
bad = sol;  bad.solinfo.RRx = sol.solinfo.x(1:ntot);
relx = rel_res(apply_c(getsol_lpivar_sop(bad,Km),F,G,nqi),apply_c(Kf,F,G,nqi),G);
fprintf('(iv) control solinfo.x: getsol(Km) vs Km(gamma, getsol(Q)) rel %.2e\n',relx);
out = struct('gam',g,'rel_b',r.rel_b,'relK',relK,'rel1',rel1,'maxK',pK,'rel_x',relx,...
    'ndec',ntot,'nx',numel(sol.solinfo.x));
end


% ========================================================================
% (iii) dispatch against legacy lpigetsol
% ========================================================================
function out = part_iii()
out = struct('nlegacy',0,'nnew',0,'nerr',0);
pvar s theta
prog = lpiprogram(s,theta,[0 1]);
[prog,P1] = poslpivar(prog,[1 1],2);                % dopvar
[prog,Q1] = lpivar(prog,[1 1;1 1],2);               % dopvar
[prog,g] = lpidecvar(prog,[2 1]);                   % dpvar
[prog,Pc] = lpivar_cdopvar(prog,struct('out',[1;1],'in',[1;1]), ...
            struct('out',{{{},{'s'}}},'in',{{{},{'s'}}}),[0 1],2);   % cdopvar in this program
N = numel(prog.decvartable);
prog.solinfo.RRx = randn(N+2,1);                    % 2 entries past decvartable, as a slack
prog.solinfo.info = struct('pinf',0,'dinf',0,'numerr',0);
Psol_leg = lpigetsol(prog,P1);                      % a fixed opvar
% 2-D program for dopvar2d / opvar2d.
pvar s1 s2 t1 t2
prog2 = lpiprogram([s1;s2],[t1;t2],[0 1;0 1]);
[prog2,P2] = lpivar(prog2,[1 1;1 1;1 1;1 1],1);     % dopvar2d
prog2.solinfo.RRx = randn(numel(prog2.decvartable),1);
prog2.solinfo.info = struct('pinf',0);
assert(isa(P2,'dopvar2d'),'(iii) lpivar did not return a dopvar2d')
P2sol = lpigetsol(prog2,P2);

leg = {  {'[]',[]}, {'double',[1 2;3 4]}, {'polynomial',s^2+3*s}, {'opvar',Psol_leg}, ...
         {'dpvar',g}, {'dpvar expr',2*g(1)+g(2)*s}, {'char',g.dvarname{1}}, {'cellstr',g.dvarname(:)'}, ...
         {'cell of dpvar',{g,2*g}}, {'dopvar',P1}, {'dopvar',Q1} };
for k = 1:numel(leg)
    [lbl,Xk] = leg{k}{:};
    a = lpigetsol(prog,Xk);     b = lpigetsol_sop(prog,Xk);
    c = subs_dvar_sop(Xk,prog.decvartable,prog.solinfo.RRx(1:N));
    assert(bit_equal(a,b),'(iii) %s: lpigetsol_sop differs from lpigetsol',lbl)
    assert(bit_equal(a,c),'(iii) %s: subs_dvar_sop differs from lpigetsol',lbl)
    out.nlegacy = out.nlegacy+1;
end
% The comparison must see a different decision vector (non-vacuous).
alt = prog;     alt.solinfo.RRx = prog.solinfo.RRx + 1;
assert(~bit_equal(lpigetsol(prog,P1),lpigetsol_sop(alt,P1)) && ...
       ~bit_equal(lpigetsol(prog,g),lpigetsol_sop(alt,g)),'(iii) bit_equal is blind to the values')
% Same for the 2-D classes, and getsol_lpivar_sop against getsol_lpivar.
leg2 = {{'dopvar2d',P2},{'opvar2d',P2sol}};
for k = 1:numel(leg2)
    [lbl,Xk] = leg2{k}{:};
    a = lpigetsol(prog2,Xk);    b = lpigetsol_sop(prog2,Xk);
    assert(bit_equal(a,b),'(iii) %s: lpigetsol_sop differs from lpigetsol',lbl)
    out.nlegacy = out.nlegacy+1;
end
assert(bit_equal(getsol_lpivar_sop(prog,P1),getsol_lpivar(prog,P1)),'(iii) getsol_lpivar_sop(dopvar) differs')
assert(bit_equal(getsol_lpivar_sop(prog2,P2),getsol_lpivar(prog2,P2)),'(iii) getsol_lpivar_sop(dopvar2d) differs')
% Names keyed, in any order: a permuted list gives the legacy answer.
pm = randperm(N);
c = subs_dvar_sop(g,prog.decvartable(pm),prog.solinfo.RRx(pm));
assert(bit_equal(c,lpigetsol(prog,g)),'(iii) permuted names change the dpvar result')

% A name listed twice in decvartable (a Gram's (i,j),(j,i)) takes the value
% of its LOWEST row, where getequation puts its At row; sosgetsol agrees.
Zc = cellstr(Pc.Zd(:));     nm1 = Zc{1};
[~,ic] = ismember(Zc,prog.decvartable);     dv = prog.solinfo.RRx(ic);
dupA = prog;    dupA.decvartable = [prog.decvartable(:); {nm1}];   dupA.solinfo.RRx = [prog.solinfo.RRx(1:N); 99];
dupB = prog;    dupB.decvartable = [{nm1}; prog.decvartable(:)];   dupB.solinfo.RRx = [99; prog.solinfo.RRx(1:N)];
dvB = dv;       dvB(1) = 99;
assert(same_ops(getsol_lpivar_sop(dupA,Pc),subs_dvar_sop(Pc,Zc,dv),0) && ...
       same_ops(getsol_lpivar_sop(dupB,Pc),subs_dvar_sop(Pc,Zc,dvB),0),'(iii) duplicate names: not the lowest row')
assert(double(sosgetsol(dupA,dpvar({nm1})))==dv(1) && double(sosgetsol(dupB,dpvar({nm1})))==99, ...
       '(iii) sosgetsol reads another row than the lowest')

% New classes.
Pcs = lpigetsol_sop(prog,Pc);
assert(isa(Pcs,'copvar') && same_ops(Pcs,getsol_lpivar_sop(prog,Pc),0),'(iii) cdopvar dispatch')
b12 = Pc.C{2,2};
assert(isa(lpigetsol_sop(prog,b12),'sopvar'),'(iii) sdopvar dispatch')
assert(isequal(lpigetsol_sop(prog,Pcs),Pcs),'(iii) copvar not returned as is')
assert(isequal(lpigetsol_sop(prog,Pcs.C{2,2}),Pcs.C{2,2}),'(iii) sopvar not returned as is')
mix = lpigetsol_sop(prog,{Pc,g;Pcs,[]});
assert(isequal(size(mix),[2 2]) && same_ops(mix{1,1},Pcs,0) && bit_equal(mix{1,2},lpigetsol(prog,g)) ...
    && isequal(mix{2,1},Pcs) && isempty(mix{2,2}),'(iii) mixed cell')
% the substitution needs names: none of Pc's are in prog2, so nothing moves
Pnone = getsol_lpivar_sop(prog2,Pc);
assert(isa(Pnone,'cdopvar') && isequal(cellstr(Pnone.Zd(:)),cellstr(Pc.Zd(:))),'(iii) absent names must stay')
out.nnew = 6;

% Errors.
bads = {struct('a',1),true,int8(3),"coeff_1",{1,struct()}};
if exist('nopvar','class')==8,  bads{end+1} = nopvar();     end
for k = 1:numel(bads)
    try
        lpigetsol_sop(prog,bads{k});
        ok = false;
    catch e
        ok = true;
        if ~iscell(bads{k})     % a cell of legacy content is legacy lpigetsol's error
            ok = strcmp(e.identifier,'lpigetsol_sop:badClass');
        end
    end
    assert(ok,'(iii) class %s not rejected',class(bads{k}))
    out.nerr = out.nerr+1;
end
unsolved = prog;    unsolved.solinfo.info = [];
try
    lpigetsol_sop(unsolved,Pc);     ok = false;
catch
    ok = true;
end
assert(ok,'(iii) unsolved program accepted');
try
    getsol_lpivar_sop(struct('decvartable',{{'a'}}),Pc);    ok = false;
catch e
    ok = strcmp(e.identifier,'getsol_lpivar_sop:notSolved');
end
assert(ok,'(iii) program without RRx accepted');
fprintf('(iii) %d legacy objects bit-identical to lpigetsol, %d new-class checks, %d rejections\n',...
    out.nlegacy,out.nnew,out.nerr);
end


% ========================================================================
% shared helpers
% ========================================================================
function tf = same_ops(A,B,tol)
% Same class, same block grid, same metadata and parameters within tol
% (absolute, max entry); tol = 0 demands bit identity.
tf = strcmp(class(A),class(B));
if ~tf,     return,     end
if isa(A,'copvar') || isa(A,'cdopvar')
    tf = isequal(size(A.C),size(B.C));
    for k = 1:numel(A.C)
        if ~tf,     return,     end
        if isempty(A.C{k}) || isempty(B.C{k})
            tf = isempty(A.C{k}) && isempty(B.C{k});
        else
            tf = same_ops(A.C{k},B.C{k},tol);
        end
    end
    if tf && isa(A,'cdopvar'),  tf = isequal(cellstr(A.Zd(:)),cellstr(B.Zd(:)));   end
    return
end
tf = isequal(A.vars,B.vars) && isequal(A.dom,B.dom) && isequal(A.dims,B.dims) ...
    && isequal(A.ZL,B.ZL) && isequal(A.ZR,B.ZR);
if ~tf,     return,     end
if isa(A,'sopvar')
    for g = 1:numel(A.params)
        tf = tf && isequal(size(A.params{g}),size(B.params{g})) ...
            && maxabs(A.params{g}-B.params{g})<=tol;
    end
else
    for g = 1:numel(A.params.A)
        tf = tf && isequal(size(A.params.A{g}),size(B.params.A{g})) ...
            && maxabs(A.params.A{g}-B.params.A{g})<=tol ...
            && isequal(size(A.params.B{g}),size(B.params.B{g})) ...
            && maxabs(A.params.B{g}-B.params.B{g})<=tol;
    end
    tf = tf && isequal(cellstr(A.Zd(:)),cellstr(B.Zd(:)));
end
end

function tf = bit_equal(a,b)
% Bit-identical data, recursively. 'isequal' cannot be used directly: the
% polynomial overload is ELEMENTWISE (SOSTOOLS400/multipoly/@polynomial/
% isequal.m) and returns an array.
tf = strcmp(class(a),class(b)) && isequal(size(a),size(b));
if ~tf,     return,     end
if isa(a,'polynomial')
    tf = isequal(a.coefficient,b.coefficient) && isequal(a.degmat,b.degmat) ...
        && isequal(a.varname,b.varname) && isequal(a.matdim,b.matdim);
elseif isobject(a)
    pr = properties(a);
    for k = 1:numel(pr)
        tf = tf && bit_equal(a.(pr{k}),b.(pr{k}));
    end
elseif isstruct(a)
    fa = fieldnames(a);
    tf = isequal(sort(fa),sort(fieldnames(b)));
    for i = 1:numel(a)
        for k = 1:numel(fa)
            tf = tf && bit_equal(a(i).(fa{k}),b(i).(fa{k}));
        end
    end
elseif iscell(a)
    for k = 1:numel(a),     tf = tf && bit_equal(a{k},b{k});    end
else
    tf = isequal(a,b);
end
end

function m = maxabs(x)
m = 0;  if ~isempty(x),     m = full(max(abs(x(:))));   end
end

function p = derangement(n)
% A permutation of 1:n with no fixed point (n >= 2).
p = 1:n;
if n<2,     return,     end
while any(p==1:n),  p = randperm(n);    end
end

function [x,w] = gauss01(n)
% Gauss-Legendre nodes and weights on [0,1] (Golub-Welsch).
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[Vv,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*Vv(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end
