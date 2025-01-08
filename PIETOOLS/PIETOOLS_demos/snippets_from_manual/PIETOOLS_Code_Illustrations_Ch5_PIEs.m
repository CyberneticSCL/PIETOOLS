% This document illustrates how PIETOOLS can be used to convert
% ODE-PDE and DDE systems to PIEs using the function "convert", as well as
% how to construct PIE systems and take interconnection using "piess" and
% "pielft".
% We refer to Chapter 5 of the manual for more context on the codes.

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
% DJ, 12/30/2024: Initial coding;
% DJ, 01/07/2025: Add some commenting;
clear

%% 5.2 Converting a PDE to a PIE
% Equations
% PDE:     \dot(x)(t,s) = x_ss(t,s),         0<s<1
% BCs:     x(t,0) + int(x,s,[0,1]) = 0,   
%          x(t,1) + int(x,s,[0,1]) = 0.
clear stateNameGenerator
pvar s t
x = pde_var(s,[0,1]);
PDE_dyn = diff(x,t) == diff(x,s,2);
PDE_BCs = [subs(x,s,0) + int(x,s,[0,1]) == 0;
           subs(x,s,1) + int(x,s,[0,1]) == 0];
PDE = [PDE_dyn; PDE_BCs];
PIE = convert(PDE,'pie')
PIE.vars
PIE.dom


%% 5.3 Converting a DDE to a PIE
% Equations
% DDE:  \dot{x}(t) = [-1.5,  0] x(t) + int_{-1}^{0} [3, 2.25] x(t+s) ds + int_{-2}^{0} [-1,  0] x(t+s) ds    
%                    [ 0.5, -1]                     [0, 0.5 ]                          [ 0, -1]  
clear DDE
DDE.A0 = [-1.5, 0; 0.5, -1];                        % Non-delayed term  A0*x(t)
DDE.Adi{1} = [3, 2.25; 0, 0.5]; DDE.tau(1) = 1;     % Distributed delay int_{-1}^{0} Ad1*x(t+s) ds
DDE.Adi{2} = [-1, 0; 0, -1];    DDE.tau(2) = 2;     % Distributed delay int_{-2}^{0} Ad1*x(t+s) ds
PIE = convert_PIETOOLS_DDE(DDE,'pie')


%% 5.4 Converting an Input-Output system to a PIE
% Equations
% PDE:     \dot(x)(t,s) = 0.5*x_ss(t,s) + s*(2-s)*w(t),         0<s<1
% Outputs: z(t) = int_{0}^{1} x(t,s) ds,
%          y(t) = x(t,1),
% BCs:     x(t,0) = u(t),   
%          x_{s}(t,1) = 0.
clear stateNameGenerator
pvar s t
pde_var state x input w output z control u sense y
x.vars = s;     x.dom = [0,1];
PDE_dyn = diff(x,t) == 0.5*diff(x,s,2) + s*(2-s)*w;
PDE_z = z == int(x,s,[0,1]);
PDE_y = y == subs(x,s,1);
PDE_BCs = [subs(x,s,0) == u; subs(diff(x,s),s,1) == 0];
PDE = [PDE_dyn; PDE_z; PDE_y; PDE_BCs];

PIE = convert(PDE,'pie')


%% 5.5.1 Declaring a simple PIE
opvar T A;
T.I = [0,1];    A.I = [0,1];
T.R.R1 = 1;  
A.R.R0 = 2;     A.R.R1 = 10;

PIE1 = piess(T,A)

%% 5.5.2 Declaring a PIE with disturbances and outputs
opvar C;
C.I = [0,1];    s = C.var1;
C.Q1 = 1-s;

PIE2 = piess(T,A,5*s,C,0)
PIE2.B1

%% 5.5.3 Declaring a PIE with actuation and sensing
opvar T C1 C2;
T.I = [0,1];        C1.I = [0,1];    C2.I = [0,1];
s = T.var1;         theta = T.var2;
T.R.R1 = -s;        T.R.R2 = -theta;
C1.Q1 = 0.5*s^2-s;  C2.Q1 = -s;

A = 0.5;    Tu = 1;     B1 = s*(2-s);
D12 = 1;    D22 = 1;

PIE3 = piess({T,0,Tu},A,{0,B1},{C1;C2},{0,D12;0,D22});

%% 5.5.4 Taking Interconnections of PIEs
% Impose the feedback law u(t)=K*v(t)
opvar K;
K.I = [0,1];    K.Q1 = 1;
PIE3_CL1 = closedLoopPIE(PIE3,K)

% Construct a Luenberger estimator with gain L*(vhat(t)-v(t))
opvar L;        s = L.var1;
L.I = [0,1];    L.Q2 = s*(1-s);
PIE3_CL2 = closedLoopPIE(PIE3,L,'observer')

% Construct an estimator with as output the control law
opvar T A C1;
T.I = [0,1];    A.I = [0,1];    C1.I = [0,1];
s = T.var1;         theta = T.var2;
T.R.R1 = -s;        T.R.R2 = -theta;
A.R.R0 = 0.5;       A.R.R1 = -s*(1-s)*theta;    A.R.R2 = A.R.R1;
C1.Q1 = 0.5*s^2-s;
PIE3_est = piess(T,A,{[],-L},{C1;K})
% Take interconnection to impose estimator based feedback
PIE3_CL3 = pielft(PIE3,PIE3_est)