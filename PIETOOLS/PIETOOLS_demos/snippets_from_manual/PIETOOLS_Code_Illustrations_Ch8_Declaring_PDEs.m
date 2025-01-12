% This document illustrates how PDEs can be declared as 'pde_struct' 
% objects in PIETOOLS using the Command Line Input (pde_var) format.
% We refer to Section 8.2 of the manual for more context on the codes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Code Illustrations
%
% Copyright (C)2024  PIETOOLS Team
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
% DJ, 01/10/2024: Initial coding;

clear;  clc;    clear stateNameGenerator

%% %% 8.2.1 Representing State, Input, and Output Variables

% Declare the spatial variables
pvar s1 s2

% Declare the PDE variables
% Option 1: full call to pde_var
x1 = pde_var('state',1,[s1;s2],[0,3;-1,1],[2;2]);
x2 = pde_var('state',1,[s1;s2],[0,3;-1,1],[2;2]);
u1 = pde_var('control',1,[],[]);
u2 = pde_var('control',1,s2,[-1,1]);
w1 = pde_var('input',1,[],[]);
w2 = pde_var('input',1,[],[]);
y1 = pde_var('sense',1,s1,[0,3]);
y2 = pde_var('sense',1,s2,[-1,1]);
z = pde_var('output',2,[],[]);

% Option 2: minimal call to pde_var
clear stateNameGenerator
x1 = pde_var([s1;s2],[0,3;-1,1]);
x2 = pde_var([s1;s2],[0,3;-1,1],[2;2]);
u1 = pde_var('control');
u2 = pde_var('control',s2,[-1,1]);
w1 = pde_var('in');
w2 = pde_var('in');
y1 = pde_var('sense',s1,[0,3]);
y2 = pde_var('sense',s2,[-1,1]);
z = pde_var('out',2);

y1



%% %% 8.2.2 Declaring Terms

%% 8.2.2a Multiplication
trm1 = 5*x1
trm1.free{1}.term{1}.C

trm2 = (3-s1)*s2*u1
trm2.C{1,1}

%% 8.2.2b Differentiation
trm3 = diff(trm1,s1,2)
trm3.free{1}.term{1}

trm4 = diff(x2,[s1;s2])

%% 8.2.2c Substitution
trm5 = subs(trm4,[s1;s2],[3;-1])
trm5.free{1}.term{1}.loc

trm6 = subs(x2,s2,-1)
trm6.free{1}.term{1}.loc

%% 8.2.2d Integration
trm7 = 10*int(x1,[s2;s1],[-1,1;0,3])
trm7.free{1}.term{1}.I{1}
trm7.free{1}.term{1}.I{2}

pvar s1_dum
trm_alt = int((s1-s1_dum)*subs(x2,s1,s1_dum),s1_dum,[0,s1])

%% 8.2.2e Delay
pvar t
trm8 = subs(u2,t,t-0.5)

%% 8.2.2f Addition
trms9 = trm2+trm3
trms9.C{1,1}
trms9.C{1,2}

trms10 = 5*diff(x1,s2,2) + trms9

%% 8.2.2 gConcatenation
RHS_x = [x1;trms10]
RHS_y = [x1+s1*w1; x2+w2]
RHS_z = [10*int(x1,[s1;s2],[0,3;-1,1]); subs(diff(x2,[s1;s2]),[s1;s2],[3;1])]



%% %% 8.2.3 Declaring Equations

%% 8.2.3a Setting state equations
pvar t
LHS1 = diff(x1,t)
LHS2 = diff(x2,'t')
LHS2.free{1}.term{1}.tdiff

x_eq1 = diff(x1,t)==x2
x_eq2 = diff(x2,'t')==trms10

x_eqs = [x_eq1; x_eq2]

%% 8.2.3b Setting output equations
y_eqs = [y1==subs(x1,s2,1)+s1*w1;
         y2==subs(x2,s1,3)+w2]
y_eqs.y{1}.term{2}

z_eqs = z==RHS_z
z_eqs.C{3,1}
z_eqs.C{3,2}

%% 8.2.3c Setting boundary conditions
BC1 = subs(x2,s2,-1)==0
BC2 = subs(x1,s1,0)==subs(u2,t,t-0.5)


% % Full equation:
pvar t
PDE = [diff(x1,t)==x2;
       diff(x2,t)==5*(diff(x1,s1,2)+diff(x1,s2,2))+(3-s1)*s2*u1;
       y1==subs(x1,s2,1)+s1*w1;
       y2==subs(x2,s1,3)+w2;
       z ==[10*int(x1,[s1;s2],[0,3;-1,1]);
            subs(diff(x2,[s1;s2]),[s1;s2],[3;1])];
       subs(x1,s1,0)==subs(u2,t,t-0.5);
       subs(x2,s1,0)==0;
       subs([x1;x2],s2,-1)==0;
       subs(diff([x1;x2],s1),s1,3)==0;
       subs(diff([x1;x2],s2),s2,1)==0]



%% %% 8.2.4 Post-Processing of PDE Structures

%% 8.2.4a Declaring controlled inputs and observed outputs
y_eqs_u = setControl(y_eqs,w1)
z_eqs_y = setObserve(z_eqs,z)

%% 8.2.4b Initializing PDE structures
PDE_i = initialize(PDE);

%% 8.2.4c Reordering components
PDE_r = reorder_comps(PDE_i);
x3 = pde_var();
PDE_2 = initialize([PDE;diff(x3,'t')==-x3]);
PDE_r2 = reorder_comps(PDE_2);

%% 8.2.4d Combining spatial variables
pvar s1 s2 s3
phi1 = pde_var(s1,[0,1]);  phi2 = pde_var(s2,[0,2]);  phi3 = pde_var(s3,[0,3]);
PDE_alt = [diff(phi1,'t')==diff(phi1,s1);
           diff(phi2,'t')==diff(phi2,s2);
           diff(phi3,'t')==diff(phi3,s3);
           subs(phi1,s1,0)==0;  subs(phi2,s2,0)==0;  subs(phi3,s3,0)==0];
PDE_alt = initialize(PDE_alt)
PDE_alt = combine_vars(PDE_alt,[0,1])

%% 8.2.4e Expanding delays
PDE_d = expand_delays(PDE_i)

%% 8.2.4f Expanding higher-order temporal derivatives
PDE_w = [diff(x1,'t',2)==5*(diff(x1,s1,2)+diff(x1,s2,2));
         subs(x1,s1,0)==0;  subs(diff(x1,s1),s1,3)==0;
         subs(x1,s2,-1)==0; subs(diff(x1,s2),s2,1)==0];
PDE_w = initialize(PDE_w)
PDE_w2 = expand_tderivatives(PDE_w)