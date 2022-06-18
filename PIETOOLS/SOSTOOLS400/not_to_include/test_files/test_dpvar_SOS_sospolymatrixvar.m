clc; clear; 

pvar th ksi x y
pvar gam

n=[3 4];
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
[prog1,V1] = DPsospolymatrixvar(prog1,Z,n);
%[prog1,V1] = DPsospolyvar(prog1,Z,n,'symmetric');
toc

disp('creating sos as pvar using sospolyvar')
tic
[prog3,V3] = sospolymatrixvar(prog3,Z,n);
toc

disp('comparing DPsospolyvar vs sospolyvar')
tic
err=dpvar2poly(V1)-V3
toc

