%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SOPVAR_MTIMES_CTRANSPOSE checks composition and adjoint of 'sopvar'
% objects against the trusted 'opvar' and 'opvar2d' classes.
%
% All checks compare against operators built and manipulated entirely in the
% established classes, so nothing here relies on the sopvar code being right.
%
% The property under test throughout is the ZL/ZR convention of the class:
% ZL holds the monomial degrees in the OUTPUT variables P.vars.out, ZR those
% in the INPUT variables P.vars.in, and a parameter P.params{k} is therefore
% of size dims(1)*numel(ZL) by dims(2)*numel(ZR). Several routines used to
% disagree with that, which stayed invisible whenever numel(ZL)==numel(ZR)
% and produced either a crash or a silently wrong answer otherwise. The
% cases below deliberately use bases of DIFFERENT sizes, and degree sets that
% differ as sets while having the same size, since the two failure modes are
% distinct.
%
% MMP, 08/29/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(17);
tol = 1e-9;
npass = 0;

pvar s1 t1 s2 s1_dum s2_dum

%% 1D conversion round trip, with unequal monomial bases
fprintf('=== 1D: sopvar2opvar(opvar2sopvar(A)) == A ===\n');
for k = 1:6
    A = rand_opvar_1d(s1,t1);
    Ab = sopvar2opvar(opvar2sopvar(A));
    assert_zero(A-Ab,tol,sprintf('1D conversion round trip, k=%d',k));
    npass = npass+1;
end
fprintf('  passed: 6 round trips\n');

%% 1D adjoint
fprintf('=== 1D: sopvar2opvar(As'') == A'' ===\n');
for k = 1:6
    A = rand_opvar_1d(s1,t1);
    got = sopvar2opvar(opvar2sopvar(A)');
    assert_zero(A'-got,tol,sprintf('1D adjoint, k=%d',k));
    npass = npass+1;
end
fprintf('  passed: 6 adjoints\n');

%% 1D composition, including adjoints on either side
fprintf('=== 1D: composition against opvar ===\n');
for k = 1:6
    A = rand_opvar_1d(s1,t1);
    B = rand_opvar_1d(s1,t1);
    As = opvar2sopvar(A);   Bs = opvar2sopvar(B);
    assert_zero(A*B   - sopvar2opvar(As*Bs),  tol,sprintf('A*B, k=%d',k));
    assert_zero(A'*B  - sopvar2opvar(As'*Bs), tol,sprintf('A''*B, k=%d',k));
    assert_zero(A*B'  - sopvar2opvar(As*Bs'), tol,sprintf('A*B'', k=%d',k));
    npass = npass+3;
end
fprintf('  passed: 18 compositions\n');

%% 2D conversion round trip and adjoint, including asymmetric spaces
fprintf('=== 2D: round trip and adjoint against opvar2d ===\n');
dom2 = [0,1;0,1];
for k = 1:10
    dim = zeros(4,2);
    dim(randi(4),1) = randi(3);
    dim(randi(4),2) = randi(3);
    A = rand_opvar2d(dim,randi(3)-1,dom2,[s1;s2],[s1_dum;s2_dum]);
    As = opvar2d2sopvar(A);
    assert_zero(A - sopvar2opvar2d(As),tol,sprintf('2D round trip, k=%d',k));
    assert_zero(A' - sopvar2opvar2d(As'),tol,sprintf('2D adjoint, k=%d',k));
    npass = npass+2;
end
fprintf('  passed: 20 round trips and adjoints\n');

%% Directly built sopvars, sweeping the monomial degree sets
fprintf('=== 1D: directly built sopvars, all degree-set pairs ===\n');
sets = { [0], [0;1], [0;1;2], [0;2], [1;2] };
dom = [0,1];
ncase = 0;
for i = 1:numel(sets)
    for j = 1:numel(sets)
        A = direct_sopvar(sets{i},sets{j},dom);
        Ao = sopvar2opvar(A);
        % adjoint
        assert_zero(Ao' - sopvar2opvar(A'),tol,...
            sprintf('direct adjoint |ZL|=%d |ZR|=%d',numel(sets{i}),numel(sets{j})));
        % compositions
        assert_zero(Ao*Ao   - sopvar2opvar(A*A),  tol,'direct A*A');
        assert_zero(Ao'*Ao  - sopvar2opvar(A'*A), tol,'direct A''*A');
        assert_zero(Ao*Ao'  - sopvar2opvar(A*A'), tol,'direct A*A''');
        % the parameter shape must follow the documented convention
        nL = A.dims(1)*prod(cellfun(@numel,A.ZL));
        nR = A.dims(2)*prod(cellfun(@numel,A.ZR));
        if ~isequal(size(A.params{1}),[nL,nR])
            error('test_sopvar_mtimes_ctranspose: parameter shape violates the ZL/ZR convention.');
        end
        At = A';
        nLt = At.dims(1)*prod(cellfun(@numel,At.ZL));
        nRt = At.dims(2)*prod(cellfun(@numel,At.ZR));
        if ~isequal(size(At.params{1}),[nLt,nRt])
            error('test_sopvar_mtimes_ctranspose: adjoint violates the ZL/ZR convention.');
        end
        ncase = ncase+1;
        npass = npass+4;
    end
end
fprintf('  passed: %d degree-set pairs, 4 checks each\n',ncase);


%% Every method must preserve the ZL/ZR convention
% This is the invariant that all of the historical swaps violated. It is
% cheap and catches the whole bug class, including for methods that have no
% convenient opvar counterpart to compare against.
fprintf('=== every method preserves the ZL/ZR convention ===\n');
A = direct_sopvar([0;1;2],[0;1],dom);       % deliberately unequal bases
B = direct_sopvar([0;2],[0;1;2],dom);
Ao = sopvar2opvar(A);   Bo = sopvar2opvar(B);

check_convention(A,'the operand itself');
check_convention(A','ctranspose');
check_convention(A*B,'mtimes');
check_convention(A+A,'plus');
check_convention([A,A],'horzcat');
check_convention([A;A],'vertcat');
check_convention(blkdiag(A,A),'blkdiag');
check_convention(repmat(A,1,2),'repmat');
check_convention(repelem(A,1,2),'repelem');
npass = npass+9;

% and where an opvar counterpart exists, the values must agree too
assert_zero((Ao+Ao) - sopvar2opvar(A+A),tol,'plus against opvar');
assert_zero([Ao,Ao] - sopvar2opvar([A,A]),tol,'horzcat against opvar');
assert_zero([Ao;Ao] - sopvar2opvar([A;A]),tol,'vertcat against opvar');
npass = npass+3;
fprintf('  passed: 9 convention checks, 3 value checks\n');

fprintf('\ntest_sopvar_mtimes_ctranspose passed (%d checks).\n',npass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_convention(P,lbl)
% ZL must have one entry per output variable, ZR one per input variable, and
% the parameters must be sized accordingly.

if numel(P.ZL)~=numel(P.vars.out)
    error('test_sopvar_mtimes_ctranspose: %s left basis has %d entries for %d output variables.',...
          lbl,numel(P.ZL),numel(P.vars.out));
end
if numel(P.ZR)~=numel(P.vars.in)
    error('test_sopvar_mtimes_ctranspose: %s right basis has %d entries for %d input variables.',...
          lbl,numel(P.ZR),numel(P.vars.in));
end
nL = P.dims(1)*prod([cellfun(@numel,P.ZL),1]);
nR = P.dims(2)*prod([cellfun(@numel,P.ZR),1]);
for k = 1:numel(P.params)
    if isempty(P.params{k})
        continue
    end
    if ~isequal(size(P.params{k}),[nL,nR])
        error('test_sopvar_mtimes_ctranspose: %s parameter %d is %s, expected %s.',...
              lbl,k,mat2str(size(P.params{k})),mat2str([nL,nR]));
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = rand_opvar_1d(var1,var2)
% A random 1D opvar whose parameters have deliberately different degrees, so
% that the left and right monomial bases end up with different sizes.

opvar A;
A.I = [0,1];
A.var1 = var1;
A.var2 = var2;
A.R.R0 = randn + randn*var1 + randn*var1^2;
A.R.R1 = randn + randn*var1 + randn*var2;
A.R.R2 = randn + randn*var1*var2 + randn*var2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = direct_sopvar(ZLdeg,ZRdeg,dom)
% Build a scalar 1D sopvar directly, following the documented convention:
% ZL over the output variable, ZR over the input variable, and the parameter
% of size numel(ZL) by numel(ZR).

vars = struct();    vars.in = {'s1'};   vars.out = {'s1'};
dd = struct();      dd.in = dom;        dd.out = dom;
params = cell(3,1);
for k = 1:3
    params{k} = sparse(randn(numel(ZLdeg),numel(ZRdeg)));
end
P = sopvar(params,vars,{ZLdeg},{ZRdeg},dd,[1,1]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assert_zero(D,tol,lbl)
% Fail unless the operator difference D is zero to the given tolerance. The
% comparison uses the operator class's own equality, after clean_opvar drops
% negligible terms, which is the idiom the existing sopvar tests use.

if ~(clean_opvar(D,tol)==0)
    error('test_sopvar_mtimes_ctranspose: %s is not zero.',lbl);
end

end
