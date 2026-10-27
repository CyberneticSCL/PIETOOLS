%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_CANONICAL_MULTIPLIER verifies the canonical multiplier form of the
% 'sopvar' and 'sdopvar' classes, and the corrected 'ctranspose' that makes
% it invariant under taking adjoints.
%
% A parameter with alpha(k)==1 carries delta(s_k-s_k'), which identifies s_k
% with s_k'. Degrees may therefore be moved between ZL and ZR in direction k
% without changing the operator, so the stored coefficients are not
% determined by it. The canonical form pins that freedom by keeping all of
% the degree on the left. Both constructors enforce it.
%
% Checks performed:
%   1. Every producer and every operation yields canonical parameters.
%   2. 'ctranspose' really is the adjoint. This is checked against the
%      defining kernel identity
%
%          Kadj_{adj(alpha)}(u,v) = K_alpha(v,u)'
%
%      evaluated from the stored parameters directly, rather than through
%      any of the machinery under test. In a multiplier direction the kernel
%      is only meaningful on the diagonal, since it is multiplied by
%      delta(s_k-s_k'), so those directions are evaluated at u_k==v_k.
%   3. Canonicalizing a non-canonical object leaves the operator alone: it
%      only moves degrees from ZR to ZL, which delta makes equivalent.
%   4. The rewrite is announced, and is idempotent.
%   5. 'eq' recognizes a self-adjoint operator as equal to its adjoint, and
%      still separates operators that genuinely differ. The first of these
%      failed before the canonical form was enforced, since taking the
%      adjoint moved degrees across and changed the stored coefficients.
%
% MMP, 08/29/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(19);
tol = 1e-9;
npass = 0;

%% 1. Every producer and operation yields canonical parameters
fprintf('=== the invariant holds for every producer and operation ===\n');
for N = 1:2
    vn = arrayfun(@(k) sprintf('s%d',k),(1:N),'UniformOutput',false);
    dm = repmat([0,1],N,1);
    vv = struct('in',{vn},'out',{vn});
    dd = struct('in',dm,'out',dm);
    gg = struct('in',ones(1,N),'out',ones(1,N));

    A = rand_sopvar([2,2],vv,dd,gg,1);
    B = rand_sopvar([2,2],vv,dd,gg,1);
    cases = { 'rand_sopvar',A ; 'ctranspose',A' ; 'plus',A+B ; 'mtimes',A*B
              'adjoint product',A'*A ; 'uminus',-A ; 'double adjoint',(A')'
              'blkdiag',blkdiag(A,B) ; 'horzcat',[A,B] ; 'vertcat',[A;B]
              'repmat',repmat(A,1,2) ; 'repelem',repelem(A,1,2) };
    for c = 1:size(cases,1)
        assert_canonical(cases{c,2},sprintf('%dD %s',N,cases{c,1}));
        npass = npass+1;
    end
end

% and for the decision variable class
for N = 1:2
    vn = arrayfun(@(k) sprintf('s%d',k),(1:N),'UniformOutput',false);
    dm = repmat([0,1],N,1);
    prog = lpiprogram(polynomial(vn(:)),polynomial(strcat(vn(:),'_dum')),dm);
    [prog,P] = possopvar(prog,2,vn,dm,1);
    for c = { 'possopvar',P ; 'ctranspose',P' ; 'plus',P+P
              'scaled',2*P ; 'self-adjoint part',P+(-1)*P' }'
        assert_canonical(c{2},sprintf('%dD sdopvar %s',N,c{1}));
        npass = npass+1;
    end
end
fprintf('  passed: every case canonical\n');

%% 2. ctranspose satisfies the defining kernel identity
fprintf('\n=== ctranspose satisfies Kadj_{adj(a)}(u,v) == K_a(v,u)'''' ===\n');
sets = { [0], [0;1], [0;1;2], [0;2], [1;2] };
worst = 0;     ncase = 0;
for i = 1:numel(sets)
    for j = 1:numel(sets)
        A = direct_sopvar(sets{i},sets{j},[0,1],[2,3]);
        worst = max(worst,adjoint_kernel_error(A,tol));
        ncase = ncase+1;
    end
end
fprintf('  passed: 1D, %d basis pairs, worst kernel error %.2e\n',ncase,worst);
npass = npass+1;

worst2 = 0;    ncase2 = 0;
for i = 1:3
    for j = 1:3
        A = direct_sopvar_2d(sets{i},sets{j},[0,1;0,1],[2,2]);
        worst2 = max(worst2,adjoint_kernel_error(A,tol));
        ncase2 = ncase2+1;
    end
end
fprintf('  passed: 2D, %d basis pairs, worst kernel error %.2e\n',ncase2,worst2);
npass = npass+1;
if ~(worst<tol && worst2<tol)
    error('test_canonical_multiplier: ctranspose does not satisfy the adjoint kernel identity.');
end

% % % The same identity where the input and output spaces differ, so that
% % % some variables carry no multiplier at all: an input-only variable
% % % becomes output-only in the adjoint and its degrees must move across,
% % % while a common variable carrying a multiplier must not. The case with
% % % no common variable at all is included, since then the parameter array
% % % is a single cell and the adjoint is a plain transpose.
fprintf('\n=== the same, with input-only and output-only variables ===\n');
mixed = { 'in {s1,s2} out {s2,s3}', {'s1','s2'}, {'s2','s3'}
          'in {s2}    out {s2,s3}', {'s2'},      {'s2','s3'}
          'in {s1,s2} out {s2}',    {'s1','s2'}, {'s2'}
          'in {s1}    out {s3}',    {'s1'},      {'s3'} };
for c = 1:size(mixed,1)
    worstm = 0;
    for trial = 1:3
        A = mixed_sopvar(mixed{c,2},mixed{c,3});
        worstm = max(worstm,adjoint_kernel_error(A,tol));
        assert_canonical(A,sprintf('%s',mixed{c,1}));
        assert_canonical(A',sprintf('%s adjoint',mixed{c,1}));
    end
    if ~(worstm<tol)
        error('test_canonical_multiplier: adjoint identity fails for %s (%.3g).',mixed{c,1},worstm);
    end
    n3 = numel(intersect(mixed{c,2},mixed{c,3}));
    fprintf('  passed: %-22s n3=%d, worst kernel error %.2e\n',mixed{c,1},n3,worstm);
    npass = npass+1;
end

%% 3. Canonicalizing does not change the operator
fprintf('\n=== canonicalization preserves the operator ===\n');
worst3 = 0;    nfold = 0;
for i = 1:numel(sets)
    for j = 1:numel(sets)
        [raw,vars,dom,dims] = raw_params(sets{i},sets{j},[0,1],[2,3]);
        [cp,ZLc,ZRc,changed] = canonicalize_multiplier(raw,vars,{sets{i}},{sets{j}},dims);
        if ~changed
            continue
        end
        nfold = nfold+1;
        e = kernel_distance(raw,{sets{i}},{sets{j}},cp,ZLc,ZRc,vars,dims);
        worst3 = max(worst3,e);
        % and the result really is canonical, and folding it again is a no-op
        [tf,~] = is_canonical_multiplier(cp,vars,ZLc,ZRc,dims);
        if ~tf
            error('test_canonical_multiplier: canonicalization did not produce a canonical object.');
        end
        [~,~,~,again] = canonicalize_multiplier(cp,vars,ZLc,ZRc,dims);
        if again
            error('test_canonical_multiplier: canonicalization is not idempotent.');
        end
    end
end
if ~(worst3<tol)
    error('test_canonical_multiplier: canonicalization changed the operator by %.3g.',worst3);
end
fprintf('  passed: %d parameter sets folded, worst kernel error %.2e, all idempotent\n',nfold,worst3);
npass = npass+2;

%% 4. The rewrite is announced
fprintf('\n=== a rewrite raises a warning ===\n');
lastwarn('');
[raw,vars,dom,~] = raw_params([0],[0;1],[0,1],[1,1]);
w = warning('off','sopvar:noncanonicalMultiplier');
cleanup = onCleanup(@() warning(w));
Pnc = sopvar(raw,vars,{[0]},{[0;1]},dom,[1,1]);   %#ok<NASGU>
[~,wid] = lastwarn;
if ~strcmp(wid,'sopvar:noncanonicalMultiplier')
    error('test_canonical_multiplier: a non-canonical input should raise sopvar:noncanonicalMultiplier, got ''%s''.',wid);
end
lastwarn('');
Pc = sopvar(Pnc.params,vars,Pnc.ZL,Pnc.ZR,dom,[1,1]);   %#ok<NASGU>
[~,wid2] = lastwarn;
if ~isempty(wid2)
    error('test_canonical_multiplier: rebuilding a canonical object should be silent, got ''%s''.',wid2);
end
fprintf('  passed: rewrite warns, rebuild of the result is silent\n');
npass = npass+2;

% % % A rewritten object has to be a usable one, not merely a correct one.
% % % An enlarged basis is built with 'union', which returns a ROW when both
% % % of its arguments are rows and treats a 1x1 array as a row, so a
% % % single-monomial ZL used to come back row-oriented. The operator was
% % % right but 'UnionBasisMonomials' read the wrong degree layout and
% % % 'plus', 'minus' and 'eq' all threw. Exercising the object is the only
% % % thing that catches that, so do it here and over every basis pair.
fprintf('\n=== a rewritten object is usable ===\n');
for op = { 'plus',@(P) P+P ; 'minus',@(P) P-P ; 'eq',@(P) P==P
           'mtimes',@(P) P*P ; 'ctranspose',@(P) P' ; 'double adjoint',@(P) (P')' }'
    try
        op{2}(Pnc);
    catch ME
        error('test_canonical_multiplier: %s fails on a rewritten object: %s',op{1},ME.message);
    end
    npass = npass+1;
end
nrow_basis = 0;
for i = 1:numel(sets)
    for j = 1:numel(sets)
        Q = direct_sopvar(sets{i},sets{j},[0,1],[1,1]);
        for r = 1:numel(Q.ZL)
            nrow_basis = nrow_basis + (size(Q.ZL{r},2)~=1);
        end
        for p = 1:numel(Q.ZR)
            nrow_basis = nrow_basis + (size(Q.ZR{p},2)~=1);
        end
        if ~(Q==Q) || ~((Q+Q)==(2*Q))
            error('test_canonical_multiplier: a rewritten object is not usable for |ZL|=%d |ZR|=%d.',...
                  numel(sets{i}),numel(sets{j}));
        end
    end
end
if nrow_basis>0
    error('test_canonical_multiplier: %d monomial bases came back row-oriented; they must be columns.',...
          nrow_basis);
end
fprintf('  passed: 6 operations on a rewritten object, %d basis pairs all column-oriented\n',...
        numel(sets)^2);
npass = npass+1;

%% 5. eq recognizes a self-adjoint operator, and still discriminates
fprintf('\n=== eq(P,adjoint(P)) for self-adjoint P ===\n');
for dg = 0:2
    prog = lpiprogram(polynomial({'s1'}),polynomial({'s1_dum'}),[0,1]);
    [prog,P] = possopvar(prog,1,{'s1'},[0,1],dg);
    if ~(P==P')
        error('test_canonical_multiplier: a possopvar operator must equal its adjoint (deg %d).',dg);
    end
    fprintf('  passed: sdopvar 1D deg %d, %d decision variables\n',dg,numel(P.Zd));
    npass = npass+1;
end
prog = lpiprogram(polynomial({'s1';'s2'}),polynomial({'s1_dum';'s2_dum'}),[0,1;0,1]);
[prog,P2] = possopvar(prog,1,{'s1','s2'},[0,1;0,1],1);
if ~(P2==P2')
    error('test_canonical_multiplier: a 2D possopvar operator must equal its adjoint.');
end
fprintf('  passed: sdopvar 2D deg 1, %d decision variables\n',numel(P2.Zd));
npass = npass+1;

% sopvar equality, which the sparse storage used to break outright
A = direct_sopvar([0;1;2],[0;1],[0,1],[2,2]);
B = direct_sopvar([0;1;2],[0;1],[0,1],[2,2]);
if ~(A==A)
    error('test_canonical_multiplier: an operator must equal itself.');
end
if ~((A')'==A)
    error('test_canonical_multiplier: the double adjoint must equal the original.');
end
if A==B
    error('test_canonical_multiplier: two independent random operators must not compare equal.');
end
if (A+A)==A
    error('test_canonical_multiplier: eq does not discriminate.');
end
fprintf('  passed: sopvar eq is reflexive, survives a double adjoint, and separates\n');
npass = npass+4;

fprintf('\ntest_canonical_multiplier passed (%d checks).\n',npass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assert_canonical(P,lbl)
% Fail unless the stored parameters satisfy the class invariant.

[tf,info] = is_canonical_multiplier(P.params,P.vars,P.ZL,P.ZR,P.dims);
if ~tf
    error('test_canonical_multiplier: %s is not canonical. %s',lbl,info.message);
end
if ~isempty(info.unchecked)
    error('test_canonical_multiplier: %s has parameters of unexpected size (%s).',...
          lbl,mat2str(info.unchecked));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = adjoint_kernel_error(A,tol)                                    %#ok<INUSD>
% Largest violation of Kadj_{adj(alpha)}(u,v) == K_alpha(v,u)', over random
% evaluation points. A multiplier direction is evaluated on the diagonal,
% since off the diagonal the kernel is multiplied by delta and carries no
% meaning.

At = A';
vS3 = intersect(A.vars.in,A.vars.out);
n3 = numel(vS3);
[~,posL] = ismember(vS3,A.vars.out);
[~,posR] = ismember(vS3,A.vars.in);
No = numel(A.vars.out);     Ni = numel(A.vars.in);
sz_C = [3*ones(1,n3),1];

e = 0;
for k = 1:numel(A.params)
    gam = ones(1,n3);
    if n3>0
        idcs = cell(1,n3);
        [idcs{:}] = ind2sub(sz_C,k);
        gam = cell2mat(idcs);
    end
    kadj = k;
    if n3>0
        adj = gam;      adj(gam==2) = 3;    adj(gam==3) = 2;
        ac = num2cell(adj);
        kadj = sub2ind(sz_C,ac{:});
    end

    for trial = 1:6
        u = 0.2+0.6*rand(1,No);     % point in the output variables of A
        v = 0.2+0.6*rand(1,Ni);     % point in the input variables of A
        % On the diagonal in every multiplier direction.
        for t = find(gam==1)
            v(posR(t)) = u(posL(t));
        end
        KA = kern(A,k,u,v);
        % The adjoint's output variables are A's input variables, and its
        % input variables are A's output variables.
        KT = kern(At,kadj,v,u);
        e = max(e,max(max(abs(KT-KA.'))));
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = kern(P,k,u,v)
% K_alpha(u,v) = (I_m o ZL(u))' C_alpha (I_n o ZR(v)), read straight off the
% stored parameter.

ZLv = monom(P.ZL,u);        ZRv = monom(P.ZR,v);
m = P.dims(1);      n = P.dims(2);
C = full(P.params{k});
if isempty(C)
    K = zeros(m,n);
    return
end
C = reshape(C,m*numel(ZLv),n*numel(ZRv));
K = kron(eye(m),ZLv)' * C * kron(eye(n),ZRv);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = monom(Z,x)
% Z(x) = x_1^Z{1} o ... o x_N^Z{N}, first variable slowest.

z = 1;
for r = 1:numel(Z)
    z = kron(z,x(r).^Z{r}(:));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = kernel_distance(p1,ZL1,ZR1,p2,ZL2,ZR2,vars,dims)
% Largest difference between the kernels of two parameter sets on the same
% variables, evaluated on the diagonal in every multiplier direction.

vS3 = intersect(vars.in,vars.out);
n3 = numel(vS3);
[~,posL] = ismember(vS3,vars.out);
[~,posR] = ismember(vS3,vars.in);
sz_C = [3*ones(1,n3),1];
S1 = struct('ZL',{ZL1},'ZR',{ZR1},'params',{p1},'dims',dims);
S2 = struct('ZL',{ZL2},'ZR',{ZR2},'params',{p2},'dims',dims);

e = 0;
for k = 1:numel(p1)
    gam = ones(1,n3);
    if n3>0
        idcs = cell(1,n3);
        [idcs{:}] = ind2sub(sz_C,k);
        gam = cell2mat(idcs);
    end
    for trial = 1:6
        u = 0.2+0.6*rand(1,numel(vars.out));
        v = 0.2+0.6*rand(1,numel(vars.in));
        for t = find(gam==1)
            v(posR(t)) = u(posL(t));
        end
        e = max(e,max(max(abs(kern(S1,k,u,v)-kern(S2,k,u,v)))));
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params,vars,dom,dims] = raw_params(ZLdeg,ZRdeg,domain,dims)
% Dense random parameters on one spatial variable, deliberately without
% regard for the canonical form.

vars = struct();    vars.in = {'s1'};   vars.out = {'s1'};
dom = struct();     dom.in = domain;    dom.out = domain;
params = cell(3,1);
for k = 1:3
    params{k} = sparse(randn(dims(1)*numel(ZLdeg),dims(2)*numel(ZRdeg)));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = direct_sopvar(ZLdeg,ZRdeg,domain,dims)
% A random 'sopvar' on one spatial variable. The constructor puts it in
% canonical form, which is what the adjoint checks then rely on.

[params,vars,dom] = raw_params(ZLdeg,ZRdeg,domain,dims);
w = warning('off','sopvar:noncanonicalMultiplier');
cleanup = onCleanup(@() warning(w));                                        %#ok<NASGU>
P = sopvar(params,vars,{ZLdeg},{ZRdeg},dom,dims);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = mixed_sopvar(vin,vout)
% A random 'sopvar' whose input and output spaces differ, with random
% monomial bases and random matrix dimensions, so that input-only,
% output-only and common variables all occur.

vars = struct('in',{vin},'out',{vout});
dom  = struct('in',repmat([0,1],numel(vin),1),'out',repmat([0,1],numel(vout),1));
ZL = arrayfun(@(k) {sort(randsample(0:3,randi([1,3]))')},1:numel(vout));
ZR = arrayfun(@(k) {sort(randsample(0:3,randi([1,3]))')},1:numel(vin));
dims = [randi([1,3]),randi([1,3])];
n3 = numel(intersect(vin,vout));
NL = prod(cellfun(@numel,ZL));      NR = prod(cellfun(@numel,ZR));
params = cell([3*ones(1,n3),1]);
for k = 1:numel(params)
    params{k} = sparse(randn(dims(1)*NL,dims(2)*NR));
end
w = warning('off','sopvar:noncanonicalMultiplier');
cleanup = onCleanup(@() warning(w));                                        %#ok<NASGU>
P = sopvar(params,vars,ZL,ZR,dom,dims);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = direct_sopvar_2d(ZLdeg,ZRdeg,domain,dims)
% A random 'sopvar' on two spatial variables, so that a parameter can carry
% a multiplier in one direction and an integral in the other.

vars = struct();    vars.in = {'s1','s2'};      vars.out = {'s1','s2'};
dom = struct();     dom.in = domain;            dom.out = domain;
nL = numel(ZLdeg)^2;    nR = numel(ZRdeg)^2;
params = cell(3,3);
for k = 1:9
    params{k} = sparse(randn(dims(1)*nL,dims(2)*nR));
end
w = warning('off','sopvar:noncanonicalMultiplier');
cleanup = onCleanup(@() warning(w));                                        %#ok<NASGU>
P = sopvar(params,vars,{ZLdeg,ZLdeg},{ZRdeg,ZRdeg},dom,dims);

end
