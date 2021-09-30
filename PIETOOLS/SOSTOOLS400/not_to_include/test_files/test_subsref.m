clc; clear; 
pvar s theta xi nu;

% inputs for left dpvar:
dvars =2; % number of decision variables
mdim=20;ndim=30; % dimensions of the dpvar matrix [m n]
dpvars=[s theta xi nu];
varname={'s', 'theta', 'xi', 'nu'};
deg=1;


% Creating first random dpvar
% processing
% construct decision variables
c=[];
for i=0:dvars-1
    c{i+1} = ['coeff_', num2str(i)];
end
dvarname=c;

nvars=length(varname);
Z2=monomials(varname,deg);
nterms=size(Z2.degmat,1);
degmat=Z2.degmat;
matdim=[mdim ndim];

% create random coefficient matrix
C = rand((dvars+1)*mdim,nterms*ndim);

% Constructing dpvar
testdpvar1 = dpvar(C,degmat,varname',dvarname,matdim);

tmp = testdpvar1(2,2:5);