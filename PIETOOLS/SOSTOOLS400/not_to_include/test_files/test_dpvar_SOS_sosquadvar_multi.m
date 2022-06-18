clc; clear; 

pvar th ksi x y
dpvar gam
n=1;
deg=1;
g=th*(1-th);
%g=th*(1-th)*ksi*(1-ksi);

mastervartable=[th x y ksi];

%mastervartable=[th x y ksi];


disp('Error testing and timing for declaration of Quadratic Matrix Variable')
Z=monomials_nd(mastervartable,0:deg);
Z1=monomials_nd([x y],0:deg);
Z2=monomials_nd([th ksi],0:deg);
nZ=size(Z,1);

prog1 = DPsosprogram(mastervartable,[gam]);

prog2 = DPsosprogram(mastervartable,[gam]);
prog3 = DPsosprogram(mastervartable,[gam]);
prog4 = DPsosprogram(mastervartable,[gam]);

disp('creating scalar sos as dpvar using DPsossosvar')
tic
[prog1,V1] = DPsossosvar(prog1,Z);
toc

disp('creating scalar sos as dpvar using DPsosquadvar, option 2')
tic
[prog2,V2] = DPsosquadvar_multi(prog2,Z,Z,1,1,'pos');
toc
err=V1-V2;

disp('creating scalar sos as dpvar using DPsosquadvar, option 1')
tic
[prog3,V3] = DPsosquadvar_opt1(prog3,Z,Z,1,1,'pos');
toc


disp('creating matrix sos as dpvar using DPsosquadvar, option 2')
prog4 = DPsosprogram(mastervartable,[gam]);
tic
[prog4,V2] = DPsosquadvar(prog4,Z,Z,3,3,'pos');
toc

disp('creating matrix sos as dpvar using DPsosquadvar, option 1')
prog4 = DPsosprogram(mastervartable,[gam]);
tic
[prog4,V2] = DPsosquadvar_opt1(prog4,Z,Z,[3 3],'pos');
toc

disp('creating matrix kernel as dpvar using DPsosquadvar, option 2')
prog4 = DPsosprogram(mastervartable,[gam]);
tic
[prog4,V2] = DPsosquadvar_multi(prog4,Z1,Z2,3, 3,'pos');
toc



disp('creating matrix kernel as dpvar using DPsosquadvar, option 1')
prog4 = DPsosprogram(mastervartable,[gam]);
tic
[prog4,V2] = DPsosquadvar_opt1(prog4,Z1,Z2,3, 3,'pos');
toc


disp('creating matrix kernel as dpvar using DPsosquadvar, option 2')
prog4 = DPsosprogram(mastervartable,[gam]);
tic
[prog4,V2] = DPsosquadvar_multi(prog4,{Z1,Z2},{Z2,Z1},[3 2], [3 2],'pos');
toc

[prog2,V2] = DPsosquadvar_multi(prog2,{th,x},{ksi,y},[1 2],[1 2],'pos');
% 
% Z2=monomials_nd([th ksi],0:deg);
% tic
% [prog1,V2] = DPsossosvar(prog1,Z2);
% toc
% 
% disp('creating sos as dpvar using DPsosposmatr')
% tic
% [prog2,P2] = DPsosposmatr(prog2,nZ);
% V=Z'*P2*Z;
% toc
% 
% disp('creating sos as pvar using sossosvar')
% tic
% [prog3,V3] = sossosvar(prog3,Z);
% toc
% 
% disp('comparing DPsossosvar vs sossosvar')
% tic
% err=dpvar2poly(V1)-V3;
% toc
% 
% disp('Constructing SOS variable using sossosmatrixvar - may take a while')
% 
% Zsym=monomials(mastervartable,2*deg);
% tic
% [prog3,P] = sospolymatrixvar(prog3,Zsym,[n n],'symmetric');
% [prog3] = sosmatrixineq(prog3,P);
% toc
% disp('Error testing and timing for construction of Positive Matrix with variables')
% disp('testing that these positive matrices are the same')
% tic
% 
% 
% 
