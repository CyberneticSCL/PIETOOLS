clc; clear; 
pvar s theta xi nu;

rng(0);

% inputs for dpvar:
dvars = 2000; % number of decision variables
mdim=4;ndim=4; % dimensions of the dpvar matrix [m n]

dpvars=[s theta xi nu];
varname={'s', 'theta', 'xi', 'nu'};
deg=4;

% inputs for right pvar:
pdim=1; % column dimension of the pvar matrix [ndim pdim]
pvars=[s theta xi nu];
pvarname={'s', 'theta', 'xi', 'nu'};
pdeg=4;
%avar=pvarname';
% ndims(avar)==2 && iscellstr(avar) && size(avar,2) == 1 

pmatdim=[ndim pdim];
pnvars=length(pvarname);
pZ2=monomials(pvars,pdeg);
pdegmat=pZ2.degmat;
pnterms=size(pdegmat,1);
pcoefficient = rand(ndim*pdim,pnterms);
testpvar=polynomial(pcoefficient',pdegmat,pvarname',pmatdim);

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

%degmat = randi([0,1],nterms,nvars)
% create random coefficient matrix
C = rand((dvars+1)*mdim,nterms*ndim);

% Constructing dpvar
testDpvar2 = dpvar(C,degmat,varname',dvarname,matdim);



%testDpvar = dpvar(c,A,[1,1]);
%testpoly = rand(1,5)*monomials(s,0:4);
%testpoly2 = c'*A; % decision polynomial without construction of dpvar

% conversion of object tyep
testDpvar_pvar=dpvar2poly(testDpvar2);
testpvar_Dpvar = poly2dpvar(testDpvar_pvar, testDpvar2.dvarname);
testpvar_Dpvar2 = poly2dpvar(testDpvar_pvar);
testDpvar_pvar-dpvar2poly(testpvar_Dpvar)
poly2dpvar(testDpvar_pvar)-testpvar_Dpvar
%%
disp('Multiplying dpvar with polynomial');
tic;
rmult1= testDpvar2*testpvar;
%lmult = testpoly*testDpvar;
toc;

%%
disp('Multiplying dpvar with polynomial as pvar');
tic;
rmult2= testDpvar_pvar*testpvar;
%lmult = testpoly*testDpvar;
toc;

% Check results
% tmp1 = dpvar2poly(rmult1);
% e12 = tmp1-rmult2;
% max(abs(e12.coefficient(:)))

%%
% Multply out the terms but implictly keep the kron(eye(mdim),Z1)
% factor on the left, i.e. don't include this into the multiplication.
% Ultimately the SOS construction will need things in factored form
% anyway.
tmp = C*kron(eye(ndim),Z2);

disp('Multiplying using pvars in a factored form');
tic;
rmulttmp= tmp*testpvar;
%lmult = testpoly*testDpvar;
toc;

% Check results
% Z1 = [1;polynomial(dvarname(:))];
% rmult3 =  kron(eye(mdim),Z1)'*rmulttmp;
% 
% e23 = rmult2-rmult3;
% max(abs(e23.coefficient(:)))

%%
% testDpvar_pvar3 = kron(eye(mdim),Z1)'*C*kron(eye(ndim),Z2)
% testDpvar_pvar

% 
% disp('Multiplying dpvar2 with polynomial');
% tic;
% rmult= testDpvar2*testpoly;
% lmult = testpoly*testDpvar2;
% toc;


% disp('Multiplying polynomial with polynomial but has same decision variables');
% tic;
%rmult2= testpoly2*testpoly;
%lmult2 = testpoly*testpoly2;
% toc;