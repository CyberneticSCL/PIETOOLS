
pvar s

% Generate random functional operators
varname1 = {'s1';'s2'};
nvars1 = numel(varname1);
varname2 = {'s2';'s3'};
nvars2 = numel(varname2);

dom = [-1,1];

m = 3;
pdeg1 = 2;
pdeg2 = 2;

Kop1 = rand_intop([m,1],varname1,dom,pdeg1);
Kop2 = rand_intop([m,1],varname2,dom,pdeg2);

% Compute the tensor product
Kop3 = otimes(Kop1,Kop2);

% Declare random state values
x1 = rand_poly([1,nvars1],s,pdeg1+1)+1;
x2 = rand_poly([1,nvars2],s,pdeg1+1)+1;

Kx1 = apply_functional(Kop1,x1,ones(1,nvars1));
Kx2 = apply_functional(Kop2,x2,ones(1,nvars2));
Kx3 = apply_functional(Kop3,[x1,x2],ones(1,nvars1+nvars2));

err = Kx3 - Kx1.*Kx2;