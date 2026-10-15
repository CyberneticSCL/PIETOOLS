% Test_tensopvar_mtimes.m -- tests @tensopvar/mtimes.m. Cases:
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
% DJ (CC), 10/15/2026: Initial coding;

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
