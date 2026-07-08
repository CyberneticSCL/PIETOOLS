clc;clear;
rng(1)
n1 = 1;
n2 = 1;
n3 = 1;

% Define sopvar
pvar s1 s2 t1 t2;
Adom = struct('in',[0,1],'out',[0,1]);
Bdom = struct('in',[0,1],'out',[0,1]);
Adims = [n1,n2]; Bdims = [n2,n3];
Avars = struct('in',{{'s1'}},'out',{{'s1'}});
Bvars = struct('in',{{'s1'}},'out',{{'s1'}});

AZL = {[0;2;3]};
AZR = {[0;1;3]};
BZ = {[0;1;2]};

Aparams = cell(3,1); Bparams = cell(3,1);

Bparams{1,1} = 1*rand(n2*3,n3*3);
Bparams{2,1} = 1*rand(n2*3,n3*3);
Bparams{3,1} = 1*rand(n2*3,n3*3);
Aparams{1,1} = 1*rand(n1*3,n2*3); 
Aparams{2,1} = 1*rand(n1*3,n2*3); 
Aparams{3,1} = 1*rand(n1*3,n2*3); 

A = sopvar(Aparams,Avars,AZL,AZR,Adom,Adims);
B = sopvar(Bparams,Bvars,BZ,BZ,Bdom,Bdims);


% convert to opvar
Aopvar= sopvar2opvar(A);
Bopvar= sopvar2opvar(B);

% opvar product
R2 = Aopvar*Bopvar;

% sopvar product
R = A*B;
R3 = sopvar2opvar(R);

diff = R2 - R3;
diff = clean_opvar(diff);
diff.R.R0
diff.R.R1
diff.R.R2 