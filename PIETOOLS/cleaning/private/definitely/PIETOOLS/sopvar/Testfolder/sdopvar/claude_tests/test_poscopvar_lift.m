%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSCOPVAR_LIFT checks the separated form Pop = A'*M*A of
% 'poscopvar_lift' ('lift_copvar' + 'posmult_cdopvar') against 'poscopvar',
% whose quadratic form it claims to reproduce at int = weight, mult = lift
% and no joint cap, and exercises what the separation adds.
%
% (1) 1-D and 2-D, mixed R x L2 and single spaces, psatz product and face,
%     sep, include, per-direction weights: the two routes have the same
%     number of decision variables, the same Gram dimension, and the same
%     linear SPAN of generators, compared in COEFFICIENT space (row g of
%     the stacked block coefficients is generator g, keyed by block, cell,
%     entry and exponents; sampled kernels misjudged 2-D ranks through
%     conditioning). The lift also passes 'verify'.
% (2) The Positivstellensatz degree count: psatz = [0 1] at weight w is
%     poscopvar at int = w plus poscopvar with the product at int = w-1;
%     a face term stays at w.
% (3) Degree vocabulary: 'int'/'mult' aliases, a scalar, per-space and
%     per-component cells (two weight groups), 'joint' with no cap
%     accepted, a capping 'joint' and 'subset' refused.
% (4) Build time of both routes, 1-D and 2-D, printed, with the three parts
%     of the separated build.
% (5) End to end: the operator V'*V + L'*L of the proof program's rank note
%     (lower kernel 2 - s), certified at lift degree 0 and weight degree 0
%     with the lower and upper components only; the solved multiplier M is
%     read back and is PSD, and the equality residual is at solver accuracy.
%
% Initial coding MMP, 10/07/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(23);

for nm = {'poscopvar','copquadvar','poscopvar_lift','lift_copvar','posmult_cdopvar','lpi_eq_sop'}
    if numel(which(nm{1},'-all'))~=1
        error('test_poscopvar_lift: %s resolves to %d files; fix the path first.',...
              nm{1},numel(which(nm{1},'-all')));
    end
end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
npass = 0;
% Parts to run: setappdata(0,'lift_test_parts',[1 3]) selects; default all.
parts = getappdata(0,'lift_test_parts');
if isempty(parts),  parts = 1:5;    end

%% (1) span equality with poscopvar at equal degrees
if any(parts==1)
% {label, dims, spaces, dom, lift d, weight w, R weight, options (both)}
cases = {
  {'1D R x L2, d=w=1',              [1;1],  {{},{'s1'}},             [0 1],      1, 1,     0, struct()}
  {'1D R^2 x L2^2, d=1 w=2',        [2;2],  {{},{'s1'}},             [0 1],      1, 2,     0, struct()}
  {'1D R x L2, d=2 w=1',            [1;1],  {{},{'s1'}},             [0 1],      2, 1,     0, struct()}
  {'1D L2 only, d=w=2',             1,      {{'s1'}},                [0 1],      2, 2,     0, struct()}
  {'1D R x L2, shifted domain',     [1;1],  {{},{'s1'}},             [-1 2],     1, 1,     0, struct()}
  {'1D R x L2, R weight 1',         [1;1],  {{},{'s1'}},             [0 1],      1, 1,     1, struct()}
  {'1D R x L2, product psatz',      [1;1],  {{},{'s1'}},             [0 1],      1, 1,     0, struct('psatz',1,'psatz_offset',0)}
  {'1D R x L2, face 3',             [1;1],  {{},{'s1'}},             [0 1],      1, 1,     0, struct('psatz',3)}
  {'1D R x L2, sep',                [1;1],  {{},{'s1'}},             [0 1],      1, 1,     0, struct('sep',true)}
  {'1D R x L2, include {[],[1;2]}', [1;1],  {{},{'s1'}},             [0 1],      1, 1,     0, struct('include',{{[],[1;2]}})}
  {'1D L2 only, d=0 w=1',           1,      {{'s1'}},                [0 1],      0, 1,     0, struct()}
  {'2D L2[s1,s2], d=w=1',           1,      {{'s1','s2'}},           [0 1;0 2],  1, 1,     0, struct()}
  {'2D R x L2[s1] x L2[s1,s2]',     [1;1;1],{{},{'s1'},{'s1','s2'}}, [0 1;0 2],  1, 1,     0, struct()}
  {'2D L2[s1] x L2[s2], w=[1 2]',   [1;1],  {{'s1'},{'s2'}},         [0 1;0 2],  1, [1 2], 0, struct()}
  {'2D L2[s1,s2], product psatz',   1,      {{'s1','s2'}},           [0 1;0 2],  1, 1,     0, struct('psatz',1,'psatz_offset',0)}
  {'2D L2[s1,s2], face 4',          1,      {{'s1','s2'}},           [0 1;0 2],  1, 1,     0, struct('psatz',4)}
  {'2D L2[s1,s2], sep [1 0]',       1,      {{'s1','s2'}},           [0 1;0 2],  1, 1,     0, struct('sep',[true,false])}
  {'2D R x L2[s1,s2], d=2 w=1',     [1;1],  {{},{'s1','s2'}},        [0 1;0 2],  2, 1,     0, struct()}
  };
for ic = 1:numel(cases)
    [lbl,dims,spaces,dom,d,w,wR,opts] = deal(cases{ic}{:});
    vars = reshape(unique([spaces{:}]),1,[]);   nv = numel(vars);
    isR = cellfun(@isempty,spaces);
    % poscopvar at the same degrees, no joint cap
    degc = cell(1,numel(spaces));
    for k = 1:numel(spaces)
        if isR(k),  degc{k} = struct('int',wR);
        else,       degc{k} = struct('int',w,'mult',d);
        end
    end
    optc = opts;
    if isfield(optc,'psatz_offset'), optc = rmfield(optc,'psatz_offset'); end
    if isfield(optc,'psatz') && optc.psatz==1 && isfield(opts,'psatz_offset') && opts.psatz_offset==0
        % product term compared at equal degree
    end
    prog = make_prog(vars,dom);
    tc0 = tic;
    [~,Pc,Qc] = poscopvar(prog,dims,spaces,dom,degc,optc);
    tcop = toc(tc0);
    % the separated form
    degl = cell(1,numel(spaces));
    for k = 1:numel(spaces)
        if isR(k),  degl{k} = struct('weight',wR);
        else,       degl{k} = struct('lift',d,'weight',w);
        end
    end
    prog = make_prog(vars,dom);
    tl0 = tic;
    [~,Pl,Ml,Al,infl] = poscopvar_lift(prog,dims,spaces,dom,degl,opts);
    tlift = toc(tl0);
    fprintf('    built %-32s poscopvar %.2fs, lift form %.2fs (compose %.2fs)\n',...
            lbl,tcop,tlift,infl.t_compose);
    v = verify(Al);
    if ~v.true
        error('test_poscopvar_lift: the lift of case ''%s'' fails verify: %s',lbl,strjoin(v.flags,' | '));
    end
    gc = sum(arrayfun(@(i) size(Qc{i,i},1),1:size(Qc,1)));
    gl = sum(infl.mult.gram);
    if numel(Pc.Zd)~=numel(Pl.Zd) || gc~=gl
        error('test_poscopvar_lift: case ''%s'': ndec %d vs %d, Gram %d vs %d.',...
              lbl,numel(Pc.Zd),numel(Pl.Zd),gc,gl);
    end
    check_no_constant(Pc,lbl);  check_no_constant(Pl,lbl);
    [r1,r2,rb] = compare_coef(Pc,Pl,vars);
    if r1~=r2 || rb~=r1
        error('test_poscopvar_lift: case ''%s'': spans differ, rank %d / %d / joint %d.',lbl,r1,r2,rb);
    end
    fprintf('  passed: %-34s ndec %4d  Gram %3d  span dim %4d\n',lbl,numel(Pl.Zd),gl,rb);
    npass = npass+1;
end
fprintf('(1) poscopvar_lift vs poscopvar: %d cases passed.\n\n',numel(cases));
end

%% (2) Positivstellensatz degree count
vars = {'s1'};  dom = [0 1];
dims = [1;1];   spaces = {{},{'s1'}};
if any(parts==2)
prog = make_prog(vars,dom);
[~,Pl] = poscopvar_lift(prog,dims,spaces,dom,struct('lift',1,'weight',2),struct('psatz',[0 1]));
% One struct applies to every space, R included (copquadvar's reading), so
% the R weight is 2 in the plain term and 1 in the product term.
prog = make_prog(vars,dom);
[prog,P0] = poscopvar(prog,dims,spaces,dom,struct('int',2,'mult',1));
[~,P1] = poscopvar(prog,dims,spaces,dom,struct('int',1,'mult',1),struct('psatz',1));
Pc = P0+P1;
[r1,r2,rb] = compare_coef(Pc,Pl,vars);
if numel(Pc.Zd)~=numel(Pl.Zd) || r1~=r2 || rb~=r1
    error('test_poscopvar_lift: psatz [0 1] at w=2 is not poscopvar(w=2) + poscopvar(product, w=1).');
end
fprintf('  passed: psatz [0 1] at w=2 = plain at 2 + product at 1 (ndec %d, span %d)\n',numel(Pl.Zd),rb);
prog = make_prog(vars,dom);
[~,Pl] = poscopvar_lift(prog,dims,spaces,dom,struct('lift',1,'weight',2),struct('psatz',3));
prog = make_prog(vars,dom);
[~,Pc] = poscopvar(prog,dims,spaces,dom,struct('int',2,'mult',1),struct('psatz',3));
[r1,r2,rb] = compare_coef(Pc,Pl,vars);
if numel(Pc.Zd)~=numel(Pl.Zd) || r1~=r2 || rb~=r1
    error('test_poscopvar_lift: a face term should stay at the weight degree.');
end
fprintf('  passed: face term at w (ndec %d, span %d)\n',numel(Pl.Zd),rb);
npass = npass+2;
fprintf('(2) Positivstellensatz degrees passed.\n\n');
end

%% (3) degree vocabulary
if any(parts==3)
prog = make_prog(vars,dom);
[~,Pa] = poscopvar_lift(prog,dims,spaces,dom,struct('int',2,'mult',1));
prog = make_prog(vars,dom);
[~,Pb] = poscopvar_lift(prog,dims,spaces,dom,struct('weight',2,'lift',1));
[r1,r2,rb] = compare_coef(Pa,Pb,vars);
if numel(Pa.Zd)~=numel(Pb.Zd) || r1~=r2 || rb~=r1
    error('test_poscopvar_lift: the int/mult aliases differ from weight/lift.');
end
fprintf('  passed: int/mult aliases\n');
% a scalar: copquadvar's reading, R weight = the scalar
prog = make_prog(vars,dom);
[~,Pa] = poscopvar_lift(prog,dims,spaces,dom,2);
prog = make_prog(vars,dom);
[~,Pb] = poscopvar(prog,dims,spaces,dom,2);
[r1,r2,rb] = compare_coef(Pa,Pb,vars);
if numel(Pa.Zd)~=numel(Pb.Zd) || r1~=r2 || rb~=r1
    error('test_poscopvar_lift: a scalar degree differs from poscopvar''s scalar degree.');
end
fprintf('  passed: scalar degree (ndec %d)\n',numel(Pa.Zd));
% joint with no cap accepted; a cap refused; subset refused
prog = make_prog(vars,dom);
[~,Pa] = poscopvar_lift(prog,dims,spaces,dom,struct('int',2,'mult',1,'joint',3));
prog = make_prog(vars,dom);
[~,Pb] = poscopvar_lift(prog,dims,spaces,dom,struct('int',2,'mult',1));
if numel(Pa.Zd)~=numel(Pb.Zd)
    error('test_poscopvar_lift: joint = int + mult should impose nothing.');
end
msg = errmsg(@() poscopvar_lift(prog,dims,spaces,dom,struct('int',2,'mult',1,'joint',2)));
if isempty(msg) || ~contains(msg,'joint')
    error('test_poscopvar_lift: a capping joint degree should be refused.');
end
msg = errmsg(@() poscopvar_lift(prog,dims,spaces,dom,struct('int',1,'subset',ones(1,4))));
if isempty(msg) || ~contains(msg,'subset')
    error('test_poscopvar_lift: a subset cap should be refused.');
end
fprintf('  passed: joint without cap accepted, capping joint and subset refused\n');
% per-component degrees: two weight groups
degl = {struct('weight',0), {struct('lift',0,'weight',1), struct('lift',1,'weight',1), struct('lift',1,'weight',2)}};
degc = {struct('int',0),    {struct('int',1,'mult',0),    struct('int',1,'mult',1),    struct('int',2,'mult',1)}};
prog = make_prog(vars,dom);
[~,Pa,~,Aa,ia] = poscopvar_lift(prog,dims,spaces,dom,degl);
prog = make_prog(vars,dom);
[~,Pb] = poscopvar(prog,dims,spaces,dom,degc);
[r1,r2,rb] = compare_coef(Pa,Pb,vars);
if numel(Pa.Zd)~=numel(Pb.Zd) || r1~=r2 || rb~=r1 || numel(ia.dims_Y)~=3
    error('test_poscopvar_lift: per-component weights (two groups) differ from poscopvar.');
end
fprintf('  passed: per-component degrees, Y = R x L2^%d x L2^%d (ndec %d)\n',...
        ia.dims_Y(2),ia.dims_Y(3),numel(Pa.Zd));
npass = npass+4;
fprintf('(3) degree vocabulary passed.\n\n');
end

%% (4) build time, both routes
if any(parts==4)
tcases = {
  {'1D R^2 x L2^2, d=2 w=2',     [2;2],  {{},{'s1'}},                       [0 1],       2, 2}
  {'1D R^2 x L2^3, d=3 w=3',     [2;3],  {{},{'s1'}},                       [0 1],       3, 3}
  {'2D L2[s1,s2], d=1 w=1',      1,      {{'s1','s2'}},                     [0 1;0 1],   1, 1}
  {'2D L2[s1,s2], d=2 w=2',      1,      {{'s1','s2'}},                     [0 1;0 1],   2, 2}
  {'2D R x L2[s1] x L2[s2] x L2[s1,s2], d=w=1', [1;1;1;1], {{},{'s1'},{'s2'},{'s1','s2'}}, [0 1;0 1], 1, 1}
  {'3D L2[s1,s2,s3], d=w=1',     1,      {{'s1','s2','s3'}},                [0 1;0 1;0 1], 1, 1}
  };
fprintf('%-44s %6s %9s %9s   %s\n','(4) build time','ndec','poscopvar','lift form','(lift / mult / compose)');
tsel = getappdata(0,'lift_test_tcases');    % e.g. 1:5 skips the 3-D case
if isempty(tsel),   tsel = 1:numel(tcases);     end
for ic = reshape(tsel,1,[])
    [lbl,dims,spaces,dom,d,w] = deal(tcases{ic}{:});
    vars = reshape(unique([spaces{:}]),1,[]);
    isR = cellfun(@isempty,spaces);
    degc = cell(1,numel(spaces));   degl = cell(1,numel(spaces));
    for k = 1:numel(spaces)
        if isR(k),  degc{k} = struct('int',0);              degl{k} = struct('weight',0);
        else,       degc{k} = struct('int',w,'mult',d);     degl{k} = struct('lift',d,'weight',w);
        end
    end
    prog = make_prog(vars,dom);
    t0 = tic;   [~,Pc] = poscopvar(prog,dims,spaces,dom,degc);      tc = toc(t0);
    prog = make_prog(vars,dom);
    t0 = tic;   [~,Pl,~,~,il] = poscopvar_lift(prog,dims,spaces,dom,degl);   tl = toc(t0);
    if numel(Pc.Zd)~=numel(Pl.Zd)
        error('test_poscopvar_lift: timing case ''%s'': ndec differ.',lbl);
    end
    fprintf('%-44s %6d %8.2fs %8.2fs   (%.2f / %.2f / %.2f)\n',lbl,numel(Pl.Zd),tc,tl,...
            il.t_lift,il.t_mult,il.t_compose);
end
fprintf('\n');
end

%% (5) end to end: V'*V + L'*L certified, the weight read back
if any(parts==5)
vars = {'s1'};  dom = [0 1];
ZL = {[0;1]};   ZR = {[0;1]};
C2 = sparse([2 0;-1 0]);    % lower kernel 2 - s   (rows s^0,s^1; cols s_dum^0,s_dum^1)
C3 = sparse([2 -1;0 0]);    % upper kernel 2 - s_dum
params = reshape({sparse(2,2),C2,C3},[3 1 1]);
Pt = copvar({sopvar(params,struct('out',{vars},'in',{vars}),ZL,ZR,...
                    struct('out',dom,'in',dom),[1 1])});
prog = lpiprogram_sop(Pt);
[prog,Pop,Mv,Av,iv] = poscopvar_lift(prog,1,{vars},dom,struct('lift',0,'weight',0),...
                                    struct('include',{{[2;3]}}));
prog = lpi_eq_sop(prog,Pop-Pt,'symmetric');
sopts = struct('solver','mosek','simplify',false);
if isempty(which('mosekopt')),  sopts.solver = 'sedumi';   end
prog = lpisolve(prog,sopts);
rel_b = cx_resid(prog);
Msol = getsol_lpivar_sop(prog,Mv);
Km = pi_blk_kernels(Msol.C{1,1},[]);
Wm = pi_poly_grid(Km{1},{'s1','s1_dum'},[0.5 0.5]);
Wm = reshape(Wm,2,2);
Psol = getsol_lpivar_sop(prog,Pop);
pts = grid_points(dom,1);
dP = norm(sample_fixed(Psol,vars,pts)-sample_fixed(Pt,vars,pts),inf);
fprintf(['(5) V''*V + L''*L at lift 0, weight 0, components {2,3}: %d decision variables, '...
         'rel_b %.1e, |Pop(sol) - P| on the grid %.1e, W = [%.3f %.3f; %.3f %.3f], eig %.3f %.3f\n'],...
        iv.ndec,rel_b,dP,Wm(1,1),Wm(1,2),Wm(2,1),Wm(2,2),min(eig(Wm)),max(eig(Wm)));
if rel_b>1e-6 || dP>1e-5 || min(eig(Wm))<-1e-6
    error('test_poscopvar_lift: the end-to-end certificate is not accepted.');
end
npass = npass+1;
end

fprintf('\ntest_poscopvar_lift passed (%d checks, parts %s).\n',npass,mat2str(parts));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prog = make_prog(vars,dom)
% An LPI program over the registry, dummies named s_dum as the classes do.
prog = lpiprogram_sop(vars,dom);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pts = grid_points(dom,nv)
% Random sample points in [s_1..s_nv, s_1_dum..s_nv_dum], each on its
% domain, on which part 5 compares the kernels of two FIXED operators.
npts = 150;     if nv>=2,   npts = 300;     end
pts = zeros(npts,2*nv);
for d = 1:nv
    pts(:,d)    = dom(d,1) + (dom(d,2)-dom(d,1))*rand(npts,1);
    pts(:,nv+d) = dom(d,1) + (dom(d,2)-dom(d,1))*rand(npts,1);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_no_constant(Pop,lbl)
% Both routes are homogeneous in the decision variables.
for i = 1:size(Pop.C,1)
    for j = 1:size(Pop.C,2)
        B = Pop.C{i,j};
        if isempty(B) || ~isa(B,'sdopvar'),     continue,   end
        for q = 1:numel(B.params.A)
            if nnz(B.params.A{q})>0
                error('test_poscopvar_lift: case ''%s'': block (%d,%d) has a constant term.',lbl,i,j);
            end
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r1,r2,rboth] = compare_coef(P1,P2,vars)
% Ranks of the generator spans of two decision containers over the same
% spaces, and of their union, in COEFFICIENT space: row g of the stacked
% B matrices of the blocks is generator g's kernel coefficients, keyed by
% (block, cell, entries, output and input exponents in registry order).
% Exact up to the floating point of the coefficients, with no evaluation
% map in between; the span tests on sampled kernels misjudged the rank by
% a few dimensions in 2-D through conditioning. The multiplier directions
% are comparable without substitution: every block is in canonical
% multiplier form (no dummy degree there), a class invariant.
[X1,k1] = coef_rows(P1,vars);
[X2,k2] = coef_rows(P2,vars);
[keys,~,loc] = unique([k1;k2],'rows');
n1 = size(k1,1);    nk = size(keys,1);
X1 = X1*sparse(1:n1,loc(1:n1),1,n1,nk);
X2 = X2*sparse(1:size(k2,1),loc(n1+1:end),1,size(k2,1),nk);
r1 = rank_rows(X1);     r2 = rank_rows(X2);
rboth = rank_rows([X1;X2]);
end


function r = rank_rows(X)
% Numerical rank of the row space. Both axes are compressed by Gaussian
% matrices to k x k before the factorization, which keeps the rank with
% probability one while it is below k (checked: on saturation k doubles),
% and the rank is read off a pivoted QR. The SVD of the row-compressed
% k x nk dense matrix was 58% of the test (profiled 10/07/2026).
[m,n] = size(X);
k = 1024;
while true
    Xc = X;
    if m>=n                 % the larger axis first, the smaller transient
        if m>k,     Xc = randn(k,m)*Xc;     end
        if n>k,     Xc = Xc*randn(n,k);     end
    else
        if n>k,     Xc = Xc*randn(n,k);     end
        if m>k,     Xc = randn(k,m)*Xc;     end
    end
    Xc = full(Xc);
    % Relative to the largest pivot: the two routes' coefficients agree to
    % rounding, so a genuine extra direction is O(1) and rounding is
    % O(1e-15); polynomial-coefficient spans have no legitimate direction
    % anywhere near 1e-8.
    [~,R,~] = qr(Xc,0);
    d = abs(diag(R));
    r = nnz(d>1e-8*max([d;0]));
    if (m>k || n>k) && r>=k,    k = 2*k;    continue,   end     % saturated
    return
end
end


function [X,keys] = coef_rows(Pop,vars)
% Generators x keyed coefficients of a decision container; a fixed block
% contributes nothing to a decision container's generators.
nv = numel(vars);
ngen = numel(Pop.Zd);
X = sparse(ngen,0);     keys = zeros(0,5+2*nv);
for i = 1:size(Pop.C,1)
    for j = 1:size(Pop.C,2)
        B = Pop.C{i,j};
        if isempty(B) || ~isa(B,'sdopvar'),     continue,   end
        m = B.dims(1);  n = B.dims(2);
        EL = multiindex_grid(B.ZL,'first_slowest');  NL = size(EL,1);
        ER = multiindex_grid(B.ZR,'first_slowest');  NR = size(ER,1);
        [~,po] = ismember(B.vars.out,vars);     [~,pi_] = ismember(B.vars.in,vars);
        eo = zeros(NL,nv);  eo(:,po) = EL;
        ei = zeros(NR,nv);  ei(:,pi_) = ER;
        % column c of B{q} <-> (rowidx, colidx), column-major, rowidx =
        % (p-1)*NL + l, colidx = (r-1)*NR + t
        nc = m*NL*n*NR;
        [ridx,cidx] = ind2sub([m*NL,n*NR],(1:nc)');
        p = floor((ridx-1)/NL)+1;   l = ridx-(p-1)*NL;
        r = floor((cidx-1)/NR)+1;   t = cidx-(r-1)*NR;
        for q = 1:numel(B.params.B)
            Bq = B.params.B{q};
            if nnz(Bq)==0,  continue,   end
            kq = [repmat([i,j,q],nc,1),p,r,eo(l,:),ei(t,:)];
            X = [X, Bq];                                                    %#ok<AGROW>
            keys = [keys; kq];                                              %#ok<AGROW>
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = sample_fixed(P,vars,pts)
% The kernels of a fixed 'copvar' on the grid, as one column.
allv = [vars, strcat(vars,'_dum')];
v = [];
for i = 1:size(P.C,1)
    for j = 1:size(P.C,2)
        B = P.C{i,j};
        if isempty(B),  continue,   end
        Kc = pi_blk_kernels(B,[]);
        shared = intersect(B.vars.out,B.vars.in);
        for q = 1:numel(Kc)
            K = polynomial(Kc{q});
            gam = gamma_of_cell(q,numel(shared));
            for t = find(gam==1)
                dn = [shared{t},'_dum'];
                if any(strcmp(K.varname,dn))
                    K = subs(K,polynomial({dn}),polynomial(shared(t)));
                end
            end
            V = pi_poly_grid(K,allv,pts);
            v = [v; V(:)];                                                      %#ok<AGROW>
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = errmsg(f)
msg = '';
try
    f();
catch ME
    msg = ME.message;
end
end
