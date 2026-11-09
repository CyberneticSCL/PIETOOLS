% Test of the 'sopvar' adjoint in two variables against 'opvar2d': for 100
% random opvar2d objects between random spaces, the adjoint of the sopvar
% conversion is converted back with sopvar2opvar2d and compared with the
% opvar2d adjoint; errors on the first mismatch.
%
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): header added;
%                  code unchanged.
for j=1:100
% Declare the independent variables
pvar s1 s2 s1_dum s2_dum
var1 = [s1;s2];
var2 = [s1_dum;s2_dum];
% Declare the domain
dom = [0,1;0,1];     % sopvar may loose domain information, need to keep [0,1]
% Declare the degree of the parameters in the independent variables
deg = randi(5)-1;
% Declare the input and output dimensions of the operators
m = randi(4);     % row dimension
n = randi(4);     % column dimension
idx1 = randi(4);    % output space
idx2 = randi(4);    % input space
dim = zeros(4,2);
dim(idx1,1) = m;
dim(idx2,2) = n;
%dim = [0,0;0,0;0,0;2,3];

% Generate random opvar2d objects of the desired dimensions
Aop = rand_opvar2d(dim,deg,dom,var1,var2);

% Compute the adjoint
Cop = Aop';

% Convert to sopvar
Asop = opvar2d2sopvar(Aop);

% Compute the sopvar adjoint
Csop = Asop';

% Compare
Cop_alt = sopvar2opvar2d(Csop);
Cop_diff = clean_opvar(Cop - Cop_alt);

if Cop_diff==0
    %disp("--- Test successful")
else
    error("--- Test failed")
end
end