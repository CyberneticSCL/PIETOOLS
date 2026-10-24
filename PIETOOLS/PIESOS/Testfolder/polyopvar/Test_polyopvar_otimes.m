% Test_polyopvar_otimes.m -- verify @polyopvar/otimes, the tensor product
% (otype=2, default) and the elementwise product (otype=1) of distributed
% polynomials, and that A.*B calls it with otype=1. Oracle: the values of
% the product, evaluated with apply_poly at random states, are compared with
% the products of the values of the factors by polynomial arithmetic: for
% the tensor product, each factor is evaluated at its own point (a copy s1_k
% of the variable), in the order of the Kronecker product; for the
% elementwise product at the same point.

% DJ (CC), 10/22/2026: Initial coding;

% NOTES
% - The case of a pointwise product as a factor of a tensor product (a
%   nested operator) needs the product of tensor-PI operators of different
%   types in @tensopvar/otimes.

tol = 1e-8;
dom = [0,1];
x1 = polyopvar('x1',polynomial({'s1'}),dom);
x2 = polyopvar('x2',polynomial({'s1'}),dom);
rng(5);
vals = struct('x1',rand_poly([1,1],{'s1'},2,4),'x2',rand_poly([1,1],{'s1'},2,4));

% vectors of linear functions of the states: scalar, 2 x 1, 3 x 1
A = [2*x1; x2+3*x1];
B = [x2; x1; x1+x2];
a = ev(A,vals);     b = ev(B,vals);
a1 = ev(x1,vals);   a2 = ev(x2,vals);

% (1) tensor product (default and otype=2) of scalars and of vectors: Kronecker order,
% each factor at its own point
check(ev(otimes(x1,x2),vals),ren(a1,1)*ren(a2,2),tol,'scalar tensor product, default type');
check(ev(otimes(x1,x2,2),vals),ren(a1,1)*ren(a2,2),tol,'scalar tensor product, otype=2');
C = otimes(A,B);
if ~isequal(size(C),[6,1])
    error("Test_polyopvar_otimes failed: the tensor product of a 2 x 1 and a 3 x 1 vector must be 6 x 1.")
end
check(ev(C,vals),kron(ren(a,1),ren(b,2)),tol,'tensor product of vectors');
check(ev(otimes(B,A),vals),kron(ren(b,1),ren(a,2)),tol,'tensor product of vectors, other order');
% a tensor product of three factors, flat
check(ev(otimes(otimes(x1,x2),x1),vals),ren(a1,1)*ren(a2,2)*ren(a1,3),tol,'tensor product of three factors');
% a pointwise product as one factor of a tensor product, on either side (nested operator)
P = x1.*x2;
pw = a1.*a2;
check(ev(otimes(P,x1),vals),ren(pw,1)*ren(a1,2),tol,'pointwise product o state');
check(ev(otimes(x1,P),vals),ren(a1,1)*ren(pw,2),tol,'state o pointwise product');

% (2) elementwise product: the same copy, sizes must match or be scalar
check(ev(otimes(A,A,1),vals),a.*a,tol,'elementwise product, otype=1');
expect_error(@() otimes(A,B,1),'must match');
% A.*B calls otimes(A,B,1)
if ~isequal(times(A,A),otimes(A,A,1)) || ~isequal(A.*A,otimes(A,A,1)) || ~isequal(x1.*x2,otimes(x1,x2,1))
    error("Test_polyopvar_otimes failed: A.*B must be equal to otimes(A,B,1).")
end
check(ev(A.*A,vals),a.*a,tol,'A.*A');

% (3) a constant
check(ev(otimes(3,B,2),vals),3*b,tol,'scalar constant, otype=2');
check(ev(otimes(3,B,1),vals),3*b,tol,'scalar constant, otype=1');
check(ev(otimes(B,3),vals),3*b,tol,'polynomial and scalar constant');
check(ev(otimes([1;2;3],B,1),vals),[1;2;3].*b,tol,'constant vector, otype=1');
expect_error(@() otimes([1;2;3],B,2),'scalar constant');
expect_error(@() otimes(B,[1;2;3],2),'scalar constant');

% (4) the type of product is a scalar 1 or 2
for ot = {[1,2],[2,2],[1;2],0,3,1.5,[],'a'}
    expect_error(@() otimes(x1,x2,ot{1}),'otype');
end

disp('Test_polyopvar_otimes passed.')


function Y = ev(P,vals)
% Value of P at the states in vals, matched by name
xv = cellfun(@(nm) vals.(nm),P.varname(:),'UniformOutput',false);
Y = apply_poly(P,xv);
end

function p = ren(p,k)
% The k-th copy s1_k of the variable s1, as named by apply_op
p.varname(strcmp(p.varname,'s1')) = {sprintf('s1_%d',k)};
end

function check(Y,E,tol,label)
err = cleanpoly(Y-E,1e-10);
if isempty(err.C) || max(abs(err.C(:)))<=tol
    return
end
error("Test_polyopvar_otimes failed (%s): residual %.3e.",label,max(abs(err.C(:))))
end

function expect_error(f,needle)
try
    f();
catch err
    if ~contains(err.message,needle)
        error('Test_polyopvar_otimes failed: wrong error message "%s" (expected "%s").',err.message,needle)
    end
    return
end
error('Test_polyopvar_otimes failed: an error containing "%s" was expected.',needle)
end
