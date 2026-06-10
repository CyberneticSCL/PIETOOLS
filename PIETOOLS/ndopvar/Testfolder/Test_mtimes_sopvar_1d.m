clc; clear;
% Declare the independent variables
pvar s1 t1
var1 = s1;
var2 = t1;
% Declare the domain
dom = [0,1];
% Declare the degree of the parameters in the independent variables
deg1 = randi(3)-1;
deg2 = randi(3)-1;
% Declare the input and output dimensions of the operators
dim1 = [0;1];
dim2 = [0;1];
dim3 = [0;1];

% Generate random opvar objects of the desired dimensions
Aop = rand_opvar([dim1,dim2],deg1,var1,var2,dom);
Bop = rand_opvar([dim2,dim3],deg2,var1,var2,dom);

% Compute the composition
Cop = Aop*Bop;

% Convert to sopvar
Asop = opvar2sopvar(Aop);
Bsop = opvar2sopvar(Bop);

% Compute the sopvar composition
Csop = Asop*Bsop;
%%
% Compare
Cop_alt = sopvar2opvar(Csop);
Cop_diff = clean_opvar(Cop - Cop_alt);
Cop_diff.R.R0
Cop_diff.R.R1
Cop_diff.R.R2

