function [Pop_test, Pop_base, Pop] = scalar_test()
% dim, deg, dom, alpha, decs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pop_scalar, Pop_base, Pop] = scalar_test(...) returns 3 nopvar objects - the 
% randomly generated nopvar (Pop), its scaled version (Pop_base.C = alpha*Pop_base.C)
% computed using cellfun, and the other scaled version (Pop_scalar) computed using
% the scalar function.

% Version: 1.0
% 
% INPUT
% dim:      1x2 array specifying the dimensions [m,n] of the operator;
% deg:      Nx1 array specifying the maximal monomial degrees d1,...,dN;
% dom:      Nx2 array with each row dom(i,:) = [ai,bi] representing the
%           spatial interval along the ith direction on which the operator
%           is defined;
% alpha:    double for multiplying the randomly generated nopvar object;
% decs:     logic variable set to true for generating ndopvar and false 
%           for nopvar;
%
% OUTPUT
% Pop:      nopvar object;
% Pop_base: scaled nopvar object computed using baseline method;
% Pop_test: scaled nopvar object computed using scale function;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding CR - 1/15/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the size of the operators (m,n), spatial dimension (N), domain (dom),
% maximal monomial degree (deg), and multiplier (alpha)
m = 1;
n = 1;
dim = [m n];
N = 2;
dom = [zeros(N,1),ones(N,1)];
deg = 2;
alpha = 2.0;
decs = true;

% Declare N variables
[N, ~] = size(dom);
var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
vars = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
dvars = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));

% Declare if decs==true decision variables
if decs
    dvarname = {'q1', 'q2', 'q3'};
else
    dvarname = cell(0,1);
end

% Declare random nopvar / ndopvar object
Pop = rand_ndopvar(dim,deg,dom,vars,dvars,dvarname);

% create ndopvar placeholder
Pop_base = Pop;

% Perform scalar multiplication on each cell element
Pop_base.C = cellfun(@(x) alpha * x, Pop.C, 'UniformOutput', false);

% computes alpha*C{ii} using scalar.
Pop_test = scalar(alpha,Pop);
 
% Check if Pop_base = Pop_scalar
if isequal(Pop_test.C, Pop_base.C)
     disp('Test passed: Scalar multiplication is correct!');
else
     disp('Test falied: Scalar multiplication is incorrect!');
end

end
