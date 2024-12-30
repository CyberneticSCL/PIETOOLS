% This document illustrates how PIETOOLS can be used to declare ODE-PDE
% system as 'pde_struct' objects using the function 'pde_var'.
% We refer to Chapter 2 of the manual for more context on the codes.


%% 4.1.1 Defining a coupled ODE-PDE system
clear stateNameGenerator
pvar t s;
x = pde_var();      
X = pde_var(s,[0,1]);
w = pde_var('in');  
z = pde_var('out',2);
u = pde_var('control');
y = pde_var('observe');

out_eq = z==[int(X,s,[0,1]); u];
eqns = [diff(x,t)==-5*x+int(diff(X,s,1),s,[0,1])+u;
        diff(X,t)==9*X+diff(X,s,2)+s*w;
        subs(X,s,0)==0;
        subs(diff(X,s),s,1)==-x+2*w;
        y==subs(X,s,0)];
odepde = [eqns; out_eq];
odepde = initialize(odepde);
PIE = convert(odepde,'pie');


%% 4.1.2 Defining 2D PDEs
clear stateNameGenerator
pvar t s1 s2
a = 0;   b = 1;   c = -1;   d = 1;
x1 = pde_var();      
x2 = pde_var(s1,[a,b]);
x3 = pde_var(s2,[c,d]);
x4 = pde_var([s1;s2],[a,b;c,d]);
w = pde_var('in',s1,[a,b]);
z = pde_var('out');
u1 = pde_var('control');     u2 = pde_var('control');
y = pde_var('observe',2,s2,[c,d]);

odepde = [diff(x1,t)==-x1+subs(x4,[s1;s2],[b;d])+u1;
         diff(x2,t)==diff(x2,s1,2)+w;
         diff(x3,t)==diff(x3,s2,2)+s2*u2;
         diff(x4,t)==diff(x4,s1,2)+diff(x4,s2,2)+4*x4;
         y==[x3;subs(x4,s1,b)];
         z==int(x3,[s1;s2],[a,b;c,d]);
         subs(x2,s1,a)==x1; subs(diff(x2,s1),s1,b)==0;
         subs(x3,s2,c)==x1; subs(x3,s2,d)==0;
         subs(x4,s2,c)==x2; subs(x4,s2,d)==0;
         subs(x4,s1,a)==x3; subs(diff(x4,s1),s1,b)==0];
odepde = initialize(odepde);
PIE = convert(odepde);


%% 4.1.3a Transport Equation
clear stateNameGenerator
pvar t s;
X = pde_var(s,[0,2]);
u = pde_var('control');  
y = pde_var('observe');
odepde = [diff(X,t)==5*diff(X,s)+u;
          subs(X,s,0)==0; 
          y==subs(X,s,2)];
odepde = initialize(odepde);
PIE = convert(odepde,'pie');


%% 4.1.3b PDE with Delay Terms
clear stateNameGenerator
pvar s t
x = pde_var();
X = pde_var(s,[0,1]);
odepde = [diff(x,t)==-5*x; 
          diff(X,t)==10*X+diff(X,s,2)+subs(x,t,t-2);
          subs(X,s,0)==0;
          subs(X,s,1)==0;];
odepde = initialize(odepde);
PIE = convert(odepde,'pie');


%% 4.1.3c Beam Equation
clear stateNameGenerator
pvar s t
x = pde_var(4,s,[0,1]);
A0 = [0,0,0,0;0,0,-1,0;0,1,0,0;0,0,0,0];
A1 = [0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0];
B = [1,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,1;0,0,0,0,0,1,0,0];
eqns = [diff(x,t)==A0*x+A1*diff(x,s); 
        B*[subs(x,s,0); subs(x,s,1)]==0];
PDE = initialize(eqns);
PIE = convert(PDE,'pie');
