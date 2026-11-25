% This document illustrates how PIETOOLS can be used to convert
% ODE-PDE and DDE systems to PIEs using the function "convert", as well as
% how to construct PIE systems and take interconnection using "piess" and
% "pielft".
% The blocks convert a heat equation with integral BCs and a wave equation
% with the same BCs, a DDE with two distributed delays, and a heat
% equation with boundary control and outputs to PIEs; declare PIEs from
% 'copvar' operators (piess), each extracted (get_op) from the PDE terms
% that it produces when applied to PDE variables; and impose a state
% feedback, a Luenberger estimator and estimator-based feedback 
% (closedLoopPIE, pielft). Lines without a semicolon print the PIE 
% structures and their fields.

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
% DJ, 01/07/2025: Add some commenting;
% DJ, 12/11/2025: Add wave equation;
% DB, MMP, DJ, 11/23/2026: Changes since PIETOOLS 2025: the operators of
%   the PIEs of 5.5 are extracted from the PDE terms that they produce when
%   when applied to a PDE variable v in L2[0,1], the constant and
%   polynomial blocks included (B1 of PIE2; A, Tu, B1, D12, D22 of PIE3); 
%   in 5.5.3-5.5.4, T has kernel -theta for theta < s and -s for theta > s, 
%   and B1 of PIE3 acts on the disturbance w; help text revised.

% Clear the workspace
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
% Convert to a PIE; print it, its spatial variables and its domain
PIE = convert(PDE,'pie')
PIE.vars
PIE.dom

% % Wave equation example:
% PDE:     \ddot(x)(t,s) = x_ss(t,s),         0<s<1
PDE_dyn = diff(x,t,2) == diff(x,s,2);
PDE = [PDE_dyn; PDE_BCs];
PIE = convert(PDE,'pie')


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
% Operators are applied to a PDE variable v in L2[0,1], and extracted from
% the resulting terms with "get_op"; v(s_dum) is v evaluated at the dummy
% variable of integration s_dum
pvar s s_dum;       dom = [0,1];                                            % DB, MMP, 11/21/2026
v = pde_var(s,dom);                                                         % DJ, 11/23/2026
T = get_op(int(v(s_dum),s_dum,[0,s]));
A = get_op(2*v + 10*int(v(s_dum),s_dum,[0,s]));

% PIE T*v_t = A*v
PIE1 = piess(T,A)

%% 5.5.2 Declaring a PIE with disturbances and outputs
% C*v = int_0^1 (1-s)*v(s) ds, and B1*w = 5*s*w for a disturbance w, as a
% function of s on the domain of v (second argument of get_op)
w = pde_var();                                                              % DJ, 11/23/2026
C = get_op(int((1-s)*v,s,[0,1]));
B1 = get_op(5*s*w,v);

% PIE T*v_t = A*v + B1*w,  z = C*v + 0*w; print it and its field B1
PIE2 = piess(T,A,B1,C,0)                                                    % DB, MMP, 11/21/2026
PIE2.B1

%% 5.5.3 Declaring a PIE with actuation and sensing
% T*v = -int_0^s s_dum*v(s_dum) ds_dum - int_s^1 s*v(s_dum) ds_dum
T = get_op(-int(s_dum*v(s_dum),s_dum,[0,s]) - int(s*v(s_dum),s_dum,[s,1])); % DJ, 11/23/2026
C1 = get_op(int((0.5*s^2-s)*v,s,[0,1]));
C2 = get_op(-int(s*v,s,[0,1]));

% Tu*u = u and B1*w = s*(2-s)*w as functions of s on the domain of v, and
% D12*u = u, for a control input u
u = pde_var();                                                              % DJ, 11/23/2026
A = get_op(0.5*v);
Tu = get_op(u,v);
B1 = get_op(s*(2-s)*w,v);
D12 = get_op(u);
D22 = D12;                                                                  % DB, MMP, 11/21/2026

% PIE d/dt (T*v + Tu*u) = A*v + B1*w + 0*u,  z = C1*v + D12*u,  y = C2*v + D22*u
PIE3 = piess({T,0,Tu},A,{B1,0},{C1;C2},{0,D12;0,D22});                      % DJ, 11/21/2026

%% 5.5.4 Taking Interconnections of PIEs
% Impose the feedback law u(t)=K*v(t)
K = get_op(int(v,s,[0,1]));                                                 % DJ, 11/23/2026
PIE3_CL1 = closedLoopPIE(PIE3,K)

% Construct a Luenberger estimator with gain L*(vhat(t)-v(t))
% L*e = s*(1-s)*e for the output error e, as a function of s on the domain
% of v
e = pde_var(1);                                                             % DJ, 11/23/2026
L = get_op(s*(1-s)*e,v);
PIE3_CL2 = closedLoopPIE(PIE3,L,'observer')

% Construct an estimator with as output the control law
T = get_op(-int(s_dum*v(s_dum),s_dum,[0,s]) - int(s*v(s_dum),s_dum,[s,1])); % DJ, 11/23/2026
A = get_op(0.5*v - int(s*(1-s)*s_dum*v(s_dum),s_dum,[0,1]));
C1 = get_op(int((0.5*s^2-s)*v,s,[0,1]));
PIE3_est = piess(T,A,{[],-L},{C1;K})
% Take interconnection to impose estimator based feedback
PIE3_CL3 = pielft(PIE3,PIE3_est)