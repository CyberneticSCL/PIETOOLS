
% Declare the independent variables
pvar s1 s2 t1 t2
var1 = [s1;s2];
var2 = [t1;t2];
% Declare the domain
dom = [-0.5,1;2,3];
% Declare the degree of the parameters in the independent variables
deg1 = randi(1)-1;
deg2 = randi(3)-1;
% Declare the input and output dimensions of the operators
dim1 = [0;0;0;1];
dim2 = [0;0;1;0];
dim3 = [0;0;1;0];

% Generate random opvar2d objects of the desired dimensions
Aop = rand_opvar2d([dim1,dim2],deg1,dom,var1,var2);
Bop = rand_opvar2d([dim2,dim3],deg2,dom,var1,var2);

% Compute the composition
Cop = Aop*Bop;

% Convert to sopvar
Asop = opvar2d2sopvar(Aop);
Bsop = opvar2d2sopvar(Bop);

% Compute the sopvar composition
Csop = Asop*Bsop;

% Compare
Cop_alt = sopvar2opvar2d(Csop);
Cop_alt.var2 = [t1;t2];
Cop_diff = clean_opvar(Cop - Cop_alt);