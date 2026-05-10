clc;clear;
pvar s1 s2 t1 t2;
Adom = struct('in',[0,1],'out',[0,1;0,1]);
Bdom = struct('in',[0,1],'out',[0,1]);
Adims = [1,1]; Bdims = [1,1];
Avars = struct('in',{{'s1'}},'out',{{'s1','s2'}});
Bvars = struct('in',{{'s1'}},'out',{{'s1'}});

AZ = {[0],[1]};
BZ = {[0;1;2]};

Aparams = cell(3,1); Bparams = cell(3,1);
for i=1:3
    Bparams{i,1} = rand(3,3);
    Aparams{i,1} = zeros(1,1);
end
Aparams{1,1} = 1; Bparams{1,1} = [rand(3,1),zeros(3,2)];

A = sopvar(Aparams,Avars,AZ,AZ,Adom,Adims);
B = sopvar(Bparams,Bvars,BZ,BZ,Bdom,Bdims);
%%
A*B