% Test of mat2sopvar. The operator is applied with 'apply_sopvar', which
% evaluates the kernel from the definition of the class, and the result is
% compared with the multiplier applied directly,
%   (Pop*x)(s) = P(s)*int x(s,v) dv,
% integrating over the input variables v that are not output variables.
% Matrices depending on a variable that is not an output variable cannot be
% multipliers, and must be refused rather than silently reduced.
%
% DJ (CC), 10/29/2026: Initial coding;
clear;
rng(1);
pvar s1 s2 s3
tol = 1e-10;
dom1 = [-1,2];

%% Multipliers between different spaces: input-only (S1), output-only (S2) and shared (S3) variables
cfgs = {struct('in',{{'s1'}},         'out',{{'s1'}}),...
        struct('in',{{'s1','s2'}},    'out',{{'s1'}}),...
        struct('in',{{'s1'}},         'out',{{'s1','s3'}}),...
        struct('in',{cell(1,0)},      'out',{{'s1','s3'}}),...
        struct('in',{{'s2','s1'}},    'out',{{'s3','s1'}}),...
        struct('in',{{'s1'}},         'out',{cell(1,0)})};
for ic=1:numel(cfgs)
    vars = cfgs{ic};
    N = numel(vars.in);     M = numel(vars.out);
    dom = struct('in',repmat(dom1,N,1),'out',repmat(dom1,M,1));
    P = rand_mat_poly(vars.out,[2,3],2);
    x = rand_poly(vars.in,3,2);
    Pop = mat2sopvar(P,vars,dom);
    check_equal(apply_sopvar(Pop,x),P*int_inputonly(x,vars,dom1),tol);
    % A matrix of type 'double' acts as a constant multiplier
    Pc = randn(2,3);
    check_equal(apply_sopvar(mat2sopvar(Pc,vars,dom),x),Pc*int_inputonly(x,vars,dom1),tol);
end

%% Matrices that depend on a variable that is not an output variable
vars = struct('in',{{'s1','s2'}},'out',{{'s1'}});
dom = struct('in',repmat(dom1,2,1),'out',dom1);
bad = {s1*s2,       's2 is an input variable only';
       s1*s3,       's3 appears in neither space';
       s1+s1*s3,    's3 in a term next to an output variable';
       1+s2,        's2 in a term next to a constant';
       [s1,s2^2],   's2 in one entry of a matrix'};
for i=1:size(bad,1)
    assert(throws(@() mat2sopvar(bad{i,1},vars,dom)),"No error when %s.",bad{i,2})
end
% A variable that cancels out is not a dependence
P = s1+s2-s2;
x = rand_poly(vars.in,1,2);
check_equal(apply_sopvar(mat2sopvar(P,vars,dom),x),s1*int_inputonly(x,vars,dom1),tol);
disp('Test_mat2sopvar passed.')



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_equal(Pix,Pexp,tol)
% CHECK_EQUAL(PIX,PEXP,TOL) checks that two polynomial vectors agree
%
% INPUTS
% - Pix:    m x 1 'polynomial', the output of the operator;
% - Pexp:   m x 1 'polynomial', the expected output;
% - tol:    scalar relative tolerance;
%
d = Pix-Pexp;
err = full(max([abs(d.C(:));0]));
scale = full(max([abs(Pexp.C(:));1]));
assert(err<=tol*scale,"Output differs from expected by %g (scale %g).",err,scale)
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function x = int_inputonly(x,vars,dom1)
% X = INT_INPUTONLY(X,VARS,DOM1) integrates X over the variables in vars.in
% that are not in vars.out, as the multiplier from L2 to a space without
% those variables acts on its input
%
% INPUTS
% - x:      n x 1 'polynomial' in the variables vars.in;
% - vars:   struct with fields 'in' and 'out', cellstr of variable names;
% - dom1:   1 x 2 'double', the domain of all variables;
%
% OUTPUTS
% - x:      n x 1 'polynomial', the integral of x over dom1 in each of the
%           variables in vars.in and not in vars.out;
%
for v = setdiff(vars.in,vars.out)
    x = int(x,polynomial(v),dom1(1),dom1(2));
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function P = rand_mat_poly(vnames,matdim,deg)
% P = RAND_MAT_POLY(VNAMES,MATDIM,DEG) returns a random polynomial matrix
%
% INPUTS
% - vnames: 1 x k cellstr, the variables of the polynomial;
% - matdim: 1 x 2 array, the row and column dimension;
% - deg:    scalar maximal degree in each variable;
%
% OUTPUTS
% - P:      matdim(1) x matdim(2) 'polynomial' in vnames with random
%           coefficients, of type 'double' if there are no variables;
%
if isempty(vnames)
    P = randn(matdim);
    return
end
pv = polynomial(reshape(vnames,[],1));
P = polynomial(randn(matdim));
for k=1:numel(vnames)
    for d=1:deg
        P = P + randn(matdim).*pv(k)^d;
    end
end
% Cross term, so the variables do not separate
P = P + randn(matdim).*prod(pv(:));
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function x = rand_poly(vnames,n,deg)
% X = RAND_POLY(VNAMES,N,DEG) returns a random polynomial vector
%
% INPUTS
% - vnames: 1 x k cellstr, the variables of the polynomial;
% - n:      scalar number of components;
% - deg:    scalar maximal degree in each variable;
%
% OUTPUTS
% - x:      n x 1 'polynomial' in vnames, with random coefficients;
%
x = polynomial(rand_mat_poly(vnames,[n,1],deg));
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function b = throws(f)
% B = THROWS(F) returns true if calling F throws an error
%
% INPUTS
% - f:      function handle without inputs;
%
% OUTPUTS
% - b:      logical, true if F() throws an error;
%
try
    f();
    b = false;
catch
    b = true;
end
end
