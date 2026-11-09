function Test_tensopvar_mtimes()
% Test_tensopvar_mtimes.m -- one local function for each case:
%   case_basic              @tensopvar/mtimes.m. Cases: (1) scalar a *
%                           multi-factor B
%   case_nested             the composition A*B of a nested 'tensopvar' A
%                           (a factor of A is itself a 'tensopvar' acting
%                           on several of the factors of the monomial, see
%                           tensopvar/plus) with an operator B
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_tensopvar_mtimes_nested into this test,
%   as cases;

cases = {@case_basic, @case_nested};
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_tensopvar_mtimes: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_tensopvar_mtimes: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- tests @tensopvar/mtimes.m. Cases:
%   (1) scalar a * multi-factor B;
%   (2) matrix * single-factor B;
%   (3) non-scalar matrix * multi-factor B must error;
%   (4) composition A*B of two tensopvars, otype(2)==2 (A) and
%       otype(1)==2 (B).
% Oracle: (1,2) apply_op(a*B,v) against a*apply_op(B,v); (4) A applied
% (apply_op) to the outputs of the factors of each term of B
% (apply_sopvar, B's terms are rank-1 Kronecker tensors), summed over the
% terms of B (linearity), against apply_op(A*B,X).
%

% NOTES
% - A non-scalar matrix multiplier is supported only for a single-factor
%   B (mtimes rescales the first factor of each term only); a scalar
%   multiplier for any number of factors;
% - Cases 2/3 use the otype(1) of the trial but have d==1 or error, so it
%   is moot. Pointwise rows (otype(1)==1) with a non-scalar multiplier
%   are unresolved (rescaling factor 1 breaks the common row count) and
%   not tested;
% - m2>=2 so that Amat is genuinely non-scalar: a 1x1 Amat is treated as
%   a scalar, which does not error;
% - Case 4: build_composable_pair makes factor j of B map onto the 'mid'
%   variable that factor j of A maps from, so that
%   A.depmat.in==B.depmat.out and A.dims(:,2)==B.dims(:,1), as mtimes
%   requires.

    ntrials = 20;
    for trial=1:ntrials
    tol = 1e-8;

    nterms = randi(3);
    N = randi(3);
    M = randi(3);
    if mod(trial,2)==1
        type = [1,2];
    else
        type = [2,2];
    end
    dnsty = 0.25;
    p = 2;
    pdegv = 2;

    %%%% Case 1: scalar multiplier, general (multi-factor) operator
    d_multi = randi(3)+1;   % d>1
    Bop_multi = build_op(d_multi,nterms,N,M,type,dnsty);
    v_multi = cell(1,d_multi);
    for j=1:d_multi
        pvars_j = polynomial(Bop_multi.ops{1,j}.vars.in(:));
        n_j = size(Bop_multi.ops{1,j},2);
        v_multi{j} = rand_poly([n_j,p],pvars_j,pdegv,50);
    end
    Bv_multi = apply_op(Bop_multi,v_multi);

    a = randn;
    C1 = a*Bop_multi;
    C1v = apply_op(C1,v_multi);
    err1 = cleanpoly(C1v - a*Bv_multi,1e-10);
    if max(max(abs(err1.C)))>tol
        error("Trial %d (type=[%d,%d]): Case 1 (a*Bop, multi-factor): residual exceeds tolerance.",trial,type(1),type(2))
    end

    %%%% Case 2: matrix multiplier, single-factor operator
    Bop_single = build_op(1,nterms,N,M,type,dnsty);
    m = size(Bop_single,1);
    v_single = cell(1,1);
    pvars_1 = polynomial(Bop_single.ops{1,1}.vars.in(:));
    n_1 = size(Bop_single.ops{1,1},2);
    v_single{1} = rand_poly([n_1,p],pvars_1,pdegv,50);
    Bv_single = apply_op(Bop_single,v_single);

    m2 = randi(2)+1;
    Amat = randn(m2,m);
    C2 = Amat*Bop_single;
    C2v = apply_op(C2,v_single);
    err2 = cleanpoly(C2v - Amat*Bv_single,1e-10);
    if max(max(abs(err2.C)))>tol
        error("Trial %d (type=[%d,%d]): Case 2 (Amat*Bop, single-factor): residual exceeds tolerance.",trial,type(1),type(2))
    end

    %%%% Case 3: matrix multiplier against a multi-factor operator must error
    threw_error = false;
    try
        C3 = Amat*Bop_multi; %#ok<NASGU>
    catch
        threw_error = true;
    end
    if ~threw_error
        error("Trial %d (type=[%d,%d]): Case 3: a non-scalar matrix multiplier against a multi-factor operator should have thrown an error.",trial,type(1),type(2))
    end

    %%%% Case 4: tensopvar*tensopvar composition (A.otype(2)==2, B.otype(1)==2)
    d4 = randi(3)+1;
    nterms_A = randi(2);
    nterms_B = randi(2);
    [Aop4,Bop4] = build_composable_pair(d4,nterms_A,nterms_B,dnsty);
    C4 = Aop4*Bop4;

    X4 = cell(1,d4);
    for j=1:d4
        pvars4_j = polynomial(Bop4.ops{1,j}.vars.in(:));
        n4_j = size(Bop4.ops{1,j},2);
        X4{j} = rand_poly([n4_j,1],pvars4_j,pdegv,50);
    end

    Y4_direct = 0;
    for i2=1:nterms_B
        Yi2 = cell(1,d4);
        for j=1:d4
            Yi2{j} = apply_sopvar(Bop4.ops{i2,j},X4{j});
        end
        Y4_direct = Y4_direct + apply_op(Aop4,Yi2);
    end
    Y4_C = apply_op(C4,X4);

    err4 = cleanpoly(Y4_C - Y4_direct,1e-10);
    if max(max(abs(err4.C)))>tol
        error("Trial %d: Case 4 (Aop*Bop, otype=[2,2] composition): residual exceeds tolerance.",trial)
    end
    end

    fprintf('Test_tensopvar_mtimes passed (%d trials, otype(1) in {1,2} alternating, incl. Case 4 composition).\n',ntrials)
end

%% --------------------------------------------------------------------------
function case_nested()
% case_nested -- tests the composition A*B of a nested
% 'tensopvar' A (a factor of A is itself a 'tensopvar' acting on several
% of the factors of the monomial, see tensopvar/plus) with an operator B
% that has one factor for each leaf factor of A (num_factors(A)), as
% subs_factor composes f.C*Rop for the substitution of some of the
% factors of a monomial.
%
% Verification is semantic and independent of mtimes: by the mixed
% product property
%   (A*B)(x1 o ... o xn) = A(B1 x1 o ... o Bn xn),
% the composition C=A*B applied to the inputs x_l of the leaf factors
% (apply_op) is compared with A applied to the outputs of the factors of
% B for these inputs (apply_sopvar, one factor at a time).
% Parts:
%   (1) the nested sum of two monomials x*m*y, B with a multiplier on the
%       functions, a gain on R, and a multiplier, with one and two terms;
%   (2) nested operators with the nesting in different factors in
%       different terms, nesting of depth two, and a nested operator with
%       several factors (multipliers from R, composed with gains on R);
%   (3) the structure of the result;
%   (4) compositions that are not supported are rejected.
%

    tol = 1e-8;
    ncases = 0;
    pvar s1 s2
    polyopvar x(s1,s2) y(s1,s2)
    Idx = x.C.ops{1,1};         Idy = y.C.ops{1,1};
    vars2 = struct('in',{{'s1','s2'}},'out',{{'s1','s2'}});
    dom2 = struct('in',[0,1;0,1],'out',[0,1;0,1]);
    varsR = struct('in',{cell(1,0)},'out',{{'s1','s2'}});
    domR = struct('in',zeros(0,2),'out',[0,1;0,1]);
    mult2 = @() mat2tensopvar(randn+randn*s1+randn*s2*s1,vars2,dom2).ops{1,1};      % multiplier on L2[0,1]^2
    gain = @() mat2tensopvar(randn).ops{1,1};                                       % gain on R
    multR = @(m) mat2tensopvar(m,varsR,domR).ops{1,1};                              % multiplier from R to L2[0,1]^2
    % (the multipliers from R depend on s1 and s2, so that nested sums of them are nested)
    m1 = 1-s1^2;        m2 = 2+s2+s1*s2;        m3 = 3-s2;

    %% (1) nested sum of x*m1*y and x*m3*y, composed with B = (multiplier) o (gain) o (multiplier)
    Ta = otimes(otimes(Idx,mat2tensopvar(m1,struct('in',{cell(1,0)},'out',{{'s1'}}),struct('in',zeros(0,2),'out',[0,1])),[1,2]),Idy,[1,2]);
    Tb = otimes(otimes(Idx,mat2tensopvar(m3,struct('in',{cell(1,0)},'out',{{'s2'}}),struct('in',zeros(0,2),'out',[0,1])),[1,2]),Idy,[1,2]);
    A = Ta+Tb;
    if ~any(cellfun('isclass',A.ops(:),'tensopvar')) || num_factors(A)~=3
        error("case_nested failed (1): setup: A should be a nested operator with 3 leaf factors.")
    end
    for KB = 1:2
        for rep = 1:3
            rng(10*KB+rep);
            B = tensopvar();        B.otype = [2,2];
            B.ops = cell(KB,3);
            for k=1:KB
                B.ops(k,:) = {mult2(),gain(),mult2()};
            end
            B = initialize(B);
            C = A*B;
            pv = polynomial({'s1','s2'}');
            X = {rand_poly([1,1],pv,2,10),randn,rand_poly([1,1],pv,2,10)};
            check_composition(A,B,C,X,tol,sprintf('(1) KB=%d, rep %d',KB,rep));
            ncases = ncases+1;
        end
    end

    %% (2) other nested operators, with multipliers from R composed with gains on R
    mR = @() multR(randn+randn*s1+randn*s2^2);
    pw = @(P,Q) pointwise(P,Q);
    kr = @(ops) kron_op(ops);
    specs = {   'term-dependent nesting',   kr({pw(mR(),mR()),mR(); mR(),pw(mR(),mR())})
                'nesting in the last factor',   kr({pw(mR(),mR()),mR(),mR(); mR(),mR(),pw(mR(),mR())})
                'depth two',               kr({pw(pw(mR(),mR()),mR()),mR(); mR(),pw(pw(mR(),mR()),mR())})
                'three terms',             kr({pw(mR(),mR()),mR(); mR(),pw(mR(),mR()); pw(mR(),mR()),mR()})
                'single nested factor',    kr({pw(mR(),mR());pw(mR(),mR())})
                'nested sum as factor',    kr({pw_terms(mR(),mR(),mR(),mR()),mR()})
                'outer pointwise product', pw(pw(mR(),mR()),mR())    };
    for k=1:size(specs,1)
        A = specs{k,2};
        nl = num_factors(A);
        for KB = 1:2
            rng(100*k+KB);
            B = tensopvar();        B.otype = [2,2];
            B.ops = cell(KB,nl);
            for t=1:KB
                for l=1:nl,     B.ops{t,l} = gain();     end
            end
            B = initialize(B);
            C = A*B;
            X = num2cell(randn(1,nl)+2);
            check_composition(A,B,C,X,tol,sprintf('(2) %s, KB=%d',specs{k,1},KB));
            ncases = ncases+1;
        end
    end

    %% (3) structure: one term for each pair of terms, the factors of A, with the leaf factors composed
    A = kr({pw(mR(),mR()),mR(); mR(),pw(mR(),mR())});
    B = tensopvar();        B.otype = [2,2];        B.ops = {gain(),gain(),gain(); gain(),gain(),gain()};   B = initialize(B);
    C = A*B;
    if ~isequal(size(C.ops),[4,2]) || num_factors(C)~=3 || ~isequal(C.otype,[A.otype(1),B.otype(2)])
        error("case_nested failed (3): size(C.ops)=%s, num_factors=%d, otype=%s.",mat2str(size(C.ops)),num_factors(C),mat2str(C.otype))
    end
    if ~isequal(cellfun('isclass',C.ops(1,:),'tensopvar'),[true,false]) || ~isequal(cellfun('isclass',C.ops(2,:),'tensopvar'),[false,true])
        error("case_nested failed (3): nested factors of the composition are in the wrong place.")
    end
    ncases = ncases+1;

    %% (4) compositions that are not supported
    A = kr({pw(mR(),mR()),mR()});               % 3 leaf factors, two factors
    Bpw = tensopvar();      Bpw.otype = [1,2];      Bpw.ops = {gain(),gain(),gain()};   Bpw = initialize(Bpw);
    expect_error(@() A*Bpw,'Kronecker');        % B cannot be split between the factors of A
    B2 = tensopvar();       B2.otype = [2,2];       B2.ops = {gain(),gain()};           B2 = initialize(B2);
    expect_error(@() A*B2,'factor');            % 2 factors for 3 leaf factors
    B4 = tensopvar();       B4.otype = [2,2];       B4.ops = {gain(),gain(),gain(),gain()};   B4 = initialize(B4);
    expect_error(@() A*B4,'factor');            % 4 factors for 3 leaf factors
    % terms with different numbers of leaf factors (the second has 4): the first sets num_factors
    Abad = kr({pw(mR(),mR()),mR(); pw(pw(mR(),mR()),mR()),mR()});
    B3 = tensopvar();       B3.otype = [2,2];       B3.ops = {gain(),gain(),gain()};    B3 = initialize(B3);
    expect_error(@() Abad*B3,'number of factors');
    % ... or fewer (the second term has 2)
    Abad2 = kr({pw(mR(),mR()),mR(); mR(),mR()});
    expect_error(@() Abad2*B3,'number of factors');
    % a pointwise product of the inputs of the factors of A (otype(2)=1): B cannot be split between the factors
    Apw_in = tensopvar();   Apw_in.otype = [2,1];   Apw_in.ops = {pw(mR(),mR()),mR()};     Apw_in = initialize(Apw_in);
    expect_error(@() Apw_in*B3,'Kronecker');
    ncases = ncases+6;

    fprintf('case_nested passed (%d cases).\n',ncases)
end

%% ---- helpers --------------------------------------------------------------

function Bop = build_op(d,nterms,N,M,type,dnsty)
% Local function must appear after all script code (MATLAB script rules).
    dims = randi(4,[d,2]);
    if type(1)==1
        dims(:,1) = dims(1,1);
    end
    if type(2)==1
        dims(:,2) = dims(1,2);
    end
    vars = struct();
    vars.in = arrayfun(@(i)sprintf('x%d',i),1:N,'UniformOutput',false);
    vars.out = arrayfun(@(i)sprintf('y%d',i),1:M,'UniformOutput',false);
    dom = struct();
    dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
    dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
    degs = struct();
    degs.in = randi(3,[1,N]);     degs.out = randi(3,[1,M]);
    depmat = struct();
    depmat.in = logical(randi([0,1],d,N));
    depmat.out = logical(randi([0,1],d,M));
    Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);
end

function [Aop,Bop] = build_composable_pair(d,nterms_A,nterms_B,dnsty)
% Builds otype=[2,2] tensopvar objects Aop,Bop with d factors each, such that Aop*Bop is well-posed:
% factor j of Bop maps some input space to a "mid" variable mid_j, and factor j of Aop maps that same
% mid_j to some output space.
    dom0 = [0,1];
    Aops = cell(nterms_A,d);        Bops = cell(nterms_B,d);
    for j=1:d
        mid_j = sprintf('mid%d',j);
        u_j = sprintf('u%d',j);         w_j = sprintf('w%d',j);
        m_u = randi(2)+1;       m_mid = randi(2)+1;      m_w = randi(2)+1;
        varsB = struct('in',{{u_j}},'out',{{mid_j}});
        for i=1:nterms_B
            Bops{i,j} = rand_sopvar([m_mid,m_u],varsB,dom0,2,dnsty);
        end
        varsA = struct('in',{{mid_j}},'out',{{w_j}});
        for i=1:nterms_A
            Aops{i,j} = rand_sopvar([m_w,m_mid],varsA,dom0,2,dnsty);
        end
    end
    Bop = tensopvar();      Bop.otype = [2,2];     Bop.ops = Bops;
    Bop = initialize(Bop);
    Aop = tensopvar();      Aop.otype = [2,2];     Aop.ops = Aops;
    Aop = initialize(Aop);
end

function check_composition(A,B,C,X,tol,label)
% Compare C=A*B applied to the leaf inputs X with A applied to the outputs of the factors of B, summed over
% the terms of B (composition is linear in B)
Yc = apply_op(C,X);
Yref = 0;
for t=1:size(B.ops,1)
    Yt = cell(1,numel(X));
    for l=1:numel(X)
        Yt{l} = apply_sopvar(B.ops{t,l},polynomial(X{l}));
    end
    Yref = Yref + apply_op(A,Yt);
end
Yc = polynomial(Yc);        Yref = polynomial(Yref);
D = Yc-Yref;
if ~isequal(size(D),size(Yref))
    error("case_nested failed %s: size %s, expected %s.",label,mat2str(size(Yc)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol*max(1,full(max(abs(nonzeros(Yref.coefficient)))))
    error("case_nested failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient)) || max(abs(nonzeros(Yref.coefficient)))<1e-6
    error("case_nested %s: reference value is zero; test would be vacuous.",label)
end
end

function N = pointwise(A,B)
N = tensopvar();        N.otype = [1,2];        N.ops = {A,B};      N = initialize(N);
end

function N = pw_terms(A1,B1,A2,B2)
% Sum of two pointwise products, a nested operator with two terms
N = tensopvar();        N.otype = [1,2];        N.ops = {A1,B1; A2,B2};     N = initialize(N);
end

function N = kron_op(ops)
N = tensopvar();        N.otype = [2,2];        N.ops = ops;        N = initialize(N);
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("case_nested failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("case_nested failed: expected an error containing '%s', but none was thrown.",str)
end
