function [Pop_test, Pop] = ctranspose_test()
% Pop_base
% dim, N, deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pop_test, Pop_base, Pop] = ctranspose_test(dim, N, deg) returns 3 nopvar objects.
% The randomly generated nopvar (Pop), its adjoint (Pop_test) and the
% baseline for comparison (Pop_base).

% Version: 1.0
% 
% INPUT
% dim:      1x2 array specifying the dimensions [m,n] of the operator;
% N:        Number of spatial variables;
% deg:      Nx1 array specifying the maximal monomial degrees d1,...,dN;
%
% OUTPUT
% Pop:      nopvar object;
% Pop_base: adjoint nopvar object computed using baseline method;
% Pop_test: scaled nopvar object computed using scale function;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding CR - 1/15/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs
m = 3;
n = 2;
dim = [m n];
N = 1;
deg = 2*ones(N,1);

% Declare domain and N variables
dom = [zeros(N,1),ones(N,1)];
var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
vars = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
dvars = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));
dvarname = cell(0,1);

% Declare random nopvar object
Pop = rand_ndopvar(dim,deg,dom,vars,dvars,dvarname);

% % create ndopvar placeholder
% Pop_base = Pop;
% 
% % Perform scalar multiplication on each cell element
% Pop_base.C = cellfun(@(x) alpha * x, Pop.C, 'UniformOutput', false);

% computes adjoint of Pop.
Pop_test = ctranspose(Pop);

end
