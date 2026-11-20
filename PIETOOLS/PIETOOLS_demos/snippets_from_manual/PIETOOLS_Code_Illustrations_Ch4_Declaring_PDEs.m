% This document illustrates how PIETOOLS can be used to declare ODE-PDE
% system as 'pde_struct' objects using the function 'pde_var'.
% Each block declares one system, initializes it and converts it to a PIE
% (the last block declares a DDE in the batch format). The variables of
% each block overwrite those of the previous one. It prints nothing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Code Illustrations
%
% Copyright (C)2026  PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, make sure to change the code in the manual as
% well, and vice versa. Document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 12/30/2024: Initial coding;
% DJ, 12/11/2025: Add wave equation example;
% DJ, 11/08/2026: Sense y=X(1) instead of X(0) in 4.1.1, as X(0)=0 by the
%                   BC made the output identically zero;
% MMP, 11/19/2026: Release comment pass (comment-only).

% Clear the workspace
clear


%% 4.1.1 Defining a coupled ODE-PDE system
% ODE state x, PDE state X on [0,1], input w, control u, outputs z, y:
%   x_t = -5x + int_0^1 X_s ds + u,     X_t = 9X + X_ss + s*w,
%   X(0) = 0,   X_s(1) = -x + 2w,   y = X(1),   z = [int_0^1 X ds; u]
clear stateNameGenerator
pvar t s;
x = pde_var();      
X = pde_var(s,[0,1]);
w = pde_var('in');  
z = pde_var('out',2);
u = pde_var('control');
y = pde_var('sense');

% Output, dynamics and BCs; initialize and convert to a PIE
out_eq = z==[int(X,s,[0,1]); u];
eqns = [diff(x,t)==-5*x+int(diff(X,s,1),s,[0,1])+u;
        diff(X,t)==9*X+diff(X,s,2)+s*w;
        subs(X,s,0)==0;
        subs(diff(X,s),s,1)==-x+2*w;
        y==subs(X,s,1)];
odepde = [eqns; out_eq];
odepde = initialize(odepde);
PIE = convert(odepde,'pie');


%% 4.1.2 Defining 2D PDEs
% ODE state x1, PDE states x2 on [a,b] (in s1), x3 on [c,d] (in s2) and
% x4 on [a,b]x[c,d]; input w(t,s1), controls u1, u2, outputs z, y(t,s2):
%   x1_t = -x1 + x4(b,d) + u1,   x2_t = x2_{s1s1} + w,   x3_t = x3_{s2s2} + s2*u2,
%   x4_t = x4_{s1s1} + x4_{s2s2} + 4*x4,   y = [x3; x4(b,s2)],
%   z = int_a^b int_c^d x3 ds2 ds1, with BCs coupling the boundary values
%   of x2, x3, x4 to x1, x2, x3
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
y = pde_var('sense',2,s2,[c,d]);

% Dynamics, outputs and BCs; initialize and convert to a PIE
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
% X_t = 5*X_s + u on [0,2], X(0) = 0, y = X(2), with the variables
% declared in command form
clear stateNameGenerator
pvar t s;
pde_var state X control u sense y;
X.vars = s;    X.dom = [0,2];
odepde = [diff(X,t)==5*diff(X,s)+u;
          subs(X,s,0)==0; 
          y==subs(X,s,2)];
odepde = initialize(odepde);
PIE = convert(odepde,'pie');


%% 4.1.3b PDE with Delay Terms
% x_t = -5x,   X_t = 10X + X_ss + x(t-2) on [0,1],   X(0) = X(1) = 0
clear stateNameGenerator
pvar s t
x = pde_var();
X = pde_var(s,[0,1]);
odepde = [diff(x,t)==-5*x; 
          diff(X,t)==10*X+diff(X,s,2)+subs(x,t,t-2);
          subs(X,s,0)==0;
          subs(X,s,1)==0];
odepde = initialize(odepde);
PIE = convert(odepde,'pie');


%% 4.1.3c Wave Equation
% x_tt = x_ss on [0,1],   x(0) = 0,   x_s(1) = 0
clear stateNameGenerator
pvar s t
x = pde_var(s,[0,1]);
eqns = [diff(x,t,2)==diff(x,s,2); 
        subs(x,s,0)==0;    subs(diff(x,s),s,1)==0];
PDE = initialize(eqns);
PIE = convert(PDE,'pie');


%% 4.1.3d Beam Equation
% Beam model in first-order form, x in L2^4[0,1]:
%   x_t = A0*x + A1*x_s,   B*[x(0); x(1)] = 0
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


%% 4.3 Batch Input Format for DDEs
% Delays tau = [1 2]:  x_t(t) = -x(t-1),   z(t) = x(t-2); INITIALIZE_PIETOOLS_DDE
% sets the remaining fields of the structure
DDE.tau = [1 2];
DDE.Ai{1} = -1;
DDE.C1i{2} = 1;
DDE = initialize_PIETOOLS_DDE(DDE);
