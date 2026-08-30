clc;clear;
pvar s1 s2 t1 t2;
Adom = struct('in',[0,1],'out',[0,1;0,1]);
Bdom = struct('in',[0,1],'out',[0,1]);
Adims = [1,1]; Bdims = [1,1];
Avars = struct('in',{{'s1'}},'out',{{'s1','s2'}});
Bvars = struct('in',{{'s1'}},'out',{{'s1'}});

AZ = {[0],[2]};
BZ = {[0;1;2]};

Aparams = cell(3,1); Bparams = cell(3,1);
for i=1:3
    Bparams{i,1} = rand(3,3);
    Aparams{i,1} = zeros(1,1);
end
Aparams{1,1} = 1; Bparams{1,1} = [rand(3,1),zeros(3,2)];

A = sopvar(Aparams,Avars,AZ(1),AZ,Adom,Adims);
B = sopvar(Bparams,Bvars,BZ,BZ,Bdom,Bdims);
%% introducing new variable diagonal test
C=A*B;


% other tests to do
% identity test  - pass
% diagonal test  - pass 
% randomized test for 1d opvar - pass
% randomized test for 2d opvar
    % sopvar to 2d conversion (and reverse)
% vector versions of the above test
% numerical integration for higher order PIs??
% other options: distributivity, associativity, 
% nopvar comparison? input/ouput spaces same, 

% 1d decompose and commute? (A\otimes B)(C\otimes D) = ??


% display(), test scripts, plus, 