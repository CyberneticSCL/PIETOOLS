clc; clear; 

pvar th ksi x y
pvar gam

deg=6;
g=th*(1-th);
%g=th*(1-th)*ksi*(1-ksi);

mastervartable=[th x y ksi];

disp('Error testing and timing for declaration of Positive Matrix Variable')
Z=monomials_nd(mastervartable,0:deg);
nZ=size(Z,1);

prog1 = sosprogram(mastervartable,[gam]);
prog2 = sosprogram(mastervartable,[gam]);
prog3 = sosprogram(mastervartable,[gam]);

disp('creating polynomial as dpvar using DPsospolyvar')
tic
[prog1,V1] = DPsospolyvar(prog1,Z);
%[prog1,V1] = DPsospolyvar(prog1,Z,'wscoeff');
toc

disp('creating sos as pvar using sospolyvar')
tic
[prog3,V3] = sospolyvar(prog3,Z);
toc

disp('comparing DPsospolyvar vs sospolyvar')
tic
err=dpvar2poly(V1)-V3;
toc

prog4 = sosprogram(mastervartable,[gam]);
Z=monomials_nd([th],0:0);
[prog4,V1] = DPsospolyvar_mat(prog4,Z,[2 3]);
dpvar2poly(V1);

prog4 = sosprogram(mastervartable,[gam]);
Z=monomials_nd([th],0:0);
[prog4,V1] = DPsospolyvar_mat(prog4,Z,[100 100],'symmetric');
dpvar2poly(V1);

 prog5 = sosprogram(mastervartable,[gam]);
 Z=monomials_nd([th],0:1,2);
 [prog5,P1] = DPsosposmatr(prog5,10);
 dpvar2poly(P1)

% prog4 = sosprogram(mastervartable,[gam]);
% Z=monomials_nd([th],0:1,2);
% [prog4,P1] = DPsosposmatr(prog4,4);
% dpvar2poly(Z'*P1*Z)
