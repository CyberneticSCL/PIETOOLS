% Test_polyopvar_native_blocks.m -- a degree-0 'polyopvar' (a known
% function) whose coefficient is a 'tensopmat' with a block of type
% 'polynomial' or 'dpvar' acting on R (see @tensopmat/tensopmat.m), as
% needed for a function with unknown coefficients, V(s)=d'*Z(s), since
% 'tensopvar' and 'sopvar' blocks cannot hold the decision variables.
% Parts:
%   (1) initialize, and the value, compared with the matrix defining the
%       block;
%   (2) sums with another degree-0 function and with a state (values
%       compared with the explicit sums);
%   (3) unary minus, scalar multiplication;
%   (4) products with a state and with another function (the block of
%       type 'polynomial' is converted to a multiplier operator for the
%       product with a state), and the error for a block of type 'dpvar'.
% The expected values are computed from the matrices defining the
% blocks, not from the class itself.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-10;
pvar s
dpvar a b
polyopvar x(s)
vars_s = struct('in',{cell(1,0)},'out',{{'s'}});    dom_S = struct('in',zeros(0,2),'out',[0,1]);
ncases = 0;

p = 1+s;        d = a*s+b;
f0 = mkdeg0(p,vars_s,dom_S);        f1 = mkdeg0(d,vars_s,dom_S);

%% (1) initialize and value
if ~f0.chkval || ~f1.chkval || ~isequal(size(f1.Z),[1,0]) || ~isequal(f1.dims.out,1) || ~isequal(f1.vars.out,{'s'})
    error("Test_polyopvar_native_blocks failed (1): degree-0 polyopvar not properly initialized.")
end
assert_equal(apply_poly(f0,{}),p,tol,'(1.1)');
assert_equal(apply_poly(f1,{}),d,tol,'(1.2)');
% vector-valued
fv = mkdeg0([1+s;a*s^2],vars_s,dom_S);
assert_equal(apply_poly(fv,{}),[1+s;a*s^2],tol,'(1.3)');
ncases = ncases+4;

%% (2) sums
f01 = f0+f0;
assert_equal(apply_poly(f01,{}),2*p,tol,'(2.1)');
f11 = f1+f1;
assert_equal(apply_poly(f11,{}),2*d,tol,'(2.2)');
f1p = f1+f0;                        % dpvar and polynomial blocks
assert_equal(apply_poly(f1p,{}),d+p,tol,'(2.3)');
xval = 2+3*s;
fx = f1+x;                          % with a state: the monomials have different types of blocks
if ~isequal(size(fx.Z,1),2)
    error("Test_polyopvar_native_blocks failed (2.4): unexpected number of monomials.")
end
assert_equal(apply_poly(fx,{xval}),d+xval,tol,'(2.4)');
fx0 = f0+3*x;
assert_equal(apply_poly(fx0,{xval}),p+3*xval,tol,'(2.5)');
ncases = ncases+5;

%% (3) unary minus, scalar multiplication
assert_equal(apply_poly(-f1,{}),-d,tol,'(3.1)');
assert_equal(apply_poly(3*f1,{}),3*d,tol,'(3.2)');
assert_equal(apply_poly(f0-f0.*1+f0,{}),p,tol,'(3.3)');
ncases = ncases+3;

%% (4) products: with a state, the block of type 'polynomial' is converted to a multiplier operator
% x*p(s) is the monomial with Z=[0,1]: the multiplier (known factor) p(s), then the state x
h = f0.*x;
if ~isequal(full(h.Z),[0,1]) || num_factors(h)~=2
    error("Test_polyopvar_native_blocks failed (4.1): Z=%s with %d factors, expected [0,1] with 2.",mat2str(full(h.Z)),num_factors(h))
end
assert_equal(apply_poly(h,{xval}),p*xval,tol,'(4.1)');
h = x.*f0;
assert_equal(apply_poly(h,{xval}),xval*p,tol,'(4.2)');
% a sum of a function and a state, times a state
g = (f0+2*x).*x;
assert_equal(apply_poly(g,{xval}),(p+2*xval)*xval,tol,'(4.3)');
% a vector-valued function, one block for each row, times the vector [x;x]: pointwise products
fv2 = mkdeg0({1+s;2*s^2},vars_s,dom_S);
h = fv2.*[x;x];
assert_equal(apply_poly(h,{xval}),[1+s;2*s^2]*xval,tol,'(4.4)');
% the product of two functions is a function: p^2
h = f0.*f0;
assert_equal(apply_poly(h,{}),p^2,tol,'(4.5)');
% a block of type 'dpvar' cannot be converted to an operator, as these cannot depend on decision variables
try
    h = f1.*x;
    error("Test_polyopvar_native_blocks failed (4.6): the product of a function with decision variables and a state should not be supported.")
catch ME
    if ~contains(ME.message,'decision variables') || contains(ME.message,'Test_polyopvar_native_blocks failed')
        error("Test_polyopvar_native_blocks failed (4.6): unexpected error message '%s'.",ME.message)
    end
end
ncases = ncases+6;

fprintf('Test_polyopvar_native_blocks passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function f = mkdeg0(P,vars,dom)
% degree-0 polyopvar with the matrix P (or cell of blocks) of type polynomial/dpvar acting on R as coefficient
f = polyopvar();        f.varname = cell(0,1);      f.Z = sparse(1,0);
if ~iscell(P),     P = {P};    end
f.C = tensopmat(P,vars,dom);
f = initialize(f);
end

function assert_equal(Y,Yref,tol,label)
if ~isequal(size(Y),size(Yref))
    error("Test_polyopvar_native_blocks failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
D = Y-Yref;
if isa(D,'double'),         c = D;
elseif isa(D,'dpvar'),      c = D.C;
else,                       c = D.coefficient;      end
err = full(max([0;abs(nonzeros(c))]));
if err>tol
    error("Test_polyopvar_native_blocks failed %s: residual %.3e.",label,err)
end
if isa(Yref,'double'),      cr = Yref;
elseif isa(Yref,'dpvar'),   cr = Yref.C;
else,                       cr = Yref.coefficient;  end
if isempty(nonzeros(cr))
    error("Test_polyopvar_native_blocks %s: reference value is zero; test would be vacuous.",label)
end
end
