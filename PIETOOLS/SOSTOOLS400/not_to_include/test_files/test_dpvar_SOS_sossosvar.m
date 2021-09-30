clc; clear; 

pvar th ksi x y
pvar gam
n=10;
deg=3;
g=th*(1-th);
%g=th*(1-th)*ksi*(1-ksi);

mastervartable=[th x y ksi];

disp('Error testing and timing for declaration of Positive Matrix Variable')
Z=monomials_nd(mastervartable,0:deg);
nZ=size(Z,1);

prog1 = sosprogram(mastervartable,[gam]);
prog2 = sosprogram(mastervartable,[gam]);
prog3 = sosprogram(mastervartable,[gam]);

disp('creating sos as dpvar using DPsossosvar')
tic
[prog1,V1] = DPsossosvar(prog1,Z);
toc
Z2=monomials_nd([th ksi],0:deg);
tic
[prog1,V2] = DPsossosvar(prog1,Z2);
toc

disp('creating sos as dpvar using DPsosposmatr')
tic
[prog2,P2] = DPsosposmatr(prog2,nZ);
V=Z'*P2*Z;
toc

disp('creating sos as pvar using sossosvar')
tic
[prog3,V3] = sossosvar(prog3,Z);
toc

disp('comparing DPsossosvar vs sossosvar')
tic
err=dpvar2poly(V1)-V3;
toc

disp('Constructing SOS variable using sossosmatrixvar - may take a while')

Zsym=monomials(mastervartable,2*deg);
tic
[prog3,P] = sospolymatrixvar(prog3,Zsym,[n n],'symmetric');
[prog3] = sosmatrixineq(prog3,P);
toc
disp('Error testing and timing for construction of Positive Matrix with variables')
disp('testing that these positive matrices are the same')
tic



