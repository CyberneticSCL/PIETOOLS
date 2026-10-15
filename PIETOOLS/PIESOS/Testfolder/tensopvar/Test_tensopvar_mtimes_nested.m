% Test_tensopvar_mtimes_nested.m -- tests the composition A*B of a nested
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
% DJ (CC), 10/15/2026: Initial coding;

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
    error("Test_tensopvar_mtimes_nested failed (1): setup: A should be a nested operator with 3 leaf factors.")
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
    error("Test_tensopvar_mtimes_nested failed (3): size(C.ops)=%s, num_factors=%d, otype=%s.",mat2str(size(C.ops)),num_factors(C),mat2str(C.otype))
end
if ~isequal(cellfun('isclass',C.ops(1,:),'tensopvar'),[true,false]) || ~isequal(cellfun('isclass',C.ops(2,:),'tensopvar'),[false,true])
    error("Test_tensopvar_mtimes_nested failed (3): nested factors of the composition are in the wrong place.")
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

fprintf('Test_tensopvar_mtimes_nested passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
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
    error("Test_tensopvar_mtimes_nested failed %s: size %s, expected %s.",label,mat2str(size(Yc)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol*max(1,full(max(abs(nonzeros(Yref.coefficient)))))
    error("Test_tensopvar_mtimes_nested failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient)) || max(abs(nonzeros(Yref.coefficient)))<1e-6
    error("Test_tensopvar_mtimes_nested %s: reference value is zero; test would be vacuous.",label)
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
        error("Test_tensopvar_mtimes_nested failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_tensopvar_mtimes_nested failed: expected an error containing '%s', but none was thrown.",str)
end
