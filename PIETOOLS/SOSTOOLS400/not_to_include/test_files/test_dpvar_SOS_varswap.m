clc; clear; 

pvar th ksi x y
pvar gam

n=[3 4];
deg=6;
g=th*(1-th);
%g=th*(1-th)*ksi*(1-ksi);
a=1;b=3; % variables to swap
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
%[prog1,V1] = DPsospolyvar(prog1,Z,n,'symmetric');
toc

disp('creating sos as pvar using sospolyvar')
tic
[prog3,V3] = sospolyvar(prog3,Z);
toc

disp('swapping variables as dpvar using DPvarswap')
tic
V1s = DPvar_swap(V1,mastervartable(a),mastervartable(b));
toc


disp('swapping variables as pvar using varswap')
tic
V3s = var_swap(V3,mastervartable(a),mastervartable(b));
toc


disp('swapping variables as pvar using subs')
tic
V3s_v = subs(V3,[mastervartable(a);mastervartable(b)],[mastervartable(b);mastervartable(a)]);
toc

disp('comparing pvar and dpvar output')
tic
err=dpvar2poly(V1s)-V3s
toc
disp('comparing pvar (subs) and dpvar output')
tic
err=dpvar2poly(V1s)-V3s_v
toc

