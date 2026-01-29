% Test script for polyopvar and tensopvar mtimes functions.

% declare two polyopvar objects.
% p1 = (A1 z1)(s1)(A2 z1)(s1)(A3 z2)(s1) + (A4 z3)(s1).
% p2 = (B1 z1)(s1)(B2 z2)(s1).

p1 = polyopvar(); p2 = polyopvar();
p1.dom=[0,1]; p2.dom=[0,1];
p1.degmat = [2,1,0;0,0,1]; p2.degmat = [1,1];
p1.varname = {'z1';'z2';'z3'}; p2.varname = {'z1';'z2'};
p1.pvarname = {'s1'}; p2.pvarname = {'s1'};

dim=[1,1]; deg=1; dom=[0,1]; var1='s1'; var2='s1_dum'; dvarname=cell(1,0);
A1 = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);

dim=[1,1]; deg=2; dom=[0,1]; var1='s1'; var2='s1_dum'; dvarname=cell(1,0);
A2 = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);

dim=[1,1]; deg=3; dom=[2,3]; var1='s1'; var2='s1_dum'; dvarname=cell(1,0);
A3 = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);

dim=[1,1]; deg=4; dom=[4,5]; var1='s1'; var2='s1_dum'; dvarname=cell(1,0);
A4 = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);

dim=[1,1]; deg=5; dom=[0,1]; var1='s1'; var2='s1_dum'; dvarname=cell(1,0);
B1 = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);

dim=[1,1]; deg=6; dom=[2,3]; var1='s1'; var2='s1_dum'; dvarname=cell(1,0);
B2 = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);

p1.C.ops = {{A1.C, A2.C, A3.C}, {A4.C}}; % nz=2; m1=m2=1; d1=3; d2=1;
p2.C.ops = {{B1.C, B2.C}}; % nz=1; m1=1; d1=2;

p1.varmat = [true;true;true]; p2.varmat = [true;true];


% Declare p3 as the expected product of polyopvars.
p3 = polyopvar();
p3.C.ops = {{A1.C,A2.C,B1.C,A3.C,B2.C}, {B1.C,B2.C,A4.C}};
p3.degmat = [3,2,0;1,1,1];
p3.varname = {'z1';'z2';'z3'};
p3.pvarname= {'s1'};
p3.dom=p1.dom;
p3.varmat=[true;true;true];


% Compute product of polyopvars.
p4 = p1*p2;

% Check polyopvars are equal - note that by assuming a common basis -
% C_nz=9 whereas C_nz=2 if calculated by hand. Problem is: how to then
% assign varname; pvarname; dom; varmat?
% varnames are same.
% pvarnames are same.
% doms are same.
% varmats are same.
% how to check if p3.C, p3.degmat are equiv. to p4.C, p4.degmat where some
% p4.C's are zeros?


