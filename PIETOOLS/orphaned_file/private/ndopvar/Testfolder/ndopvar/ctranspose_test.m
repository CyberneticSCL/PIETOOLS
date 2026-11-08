function [Pop, Pop_base, Pop_test] = ctranspose_test(dim,deg) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pop_test, Pop_base, Pop] = ctranspose_test(dim,deg) returns 3 nopvar objects.
% The randomly generated nopvar (Pop), its adjoint (Pop_test) and the
% baseline for comparison (Pop_base).

% Version: 1.0
% 
% INPUT
% dim:      1x2 array specifying the dimensions [m,n] of the operator;
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
N = 2; % must be 2 to compare with opvar2d

% Declare domain and N variables
dom = [zeros(N,1),ones(N,1)];
var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
vars = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
dvars = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));

% Declare random opvar2d object
Pop = rand_opvar2d(dim,deg,dom,vars,dvars);

% set appropriate elements of Pop.R22 to zero such that it is a 4 PI operator;
for i=1:3
    Pop.R22{1,i} = 0; Pop.R22{i,1} = 0;
end

% computes adjoint on opvar2d then convert to nopvar.
Pop_base = ctranspose(Pop);
Pop_base = dopvar2d2ndopvar(Pop_base,deg);

% Compute adjoint using nopvar ctranspose.
Pop = dopvar2d2ndopvar(Pop,deg);
Pop_test = ctranspose(Pop);

% Check if Pop_base = Pop_test
if isequal(Pop_test.C(2:3,2:3), Pop_base.C(2:3,2:3))
     disp('Test passed: ctranspose is correct!');
else
     disp('Test falied: ctranspose is incorrect!');
end

end
