% This document illustrates how PIETOOLS can be used to convert
% ODE-PDE and DDE systems to PIEs using the function "convert", as well as
% how to construct PIE systems and take interconnection using "piess" and
% "pielft".
% We refer to Chapter 5 of the manual for more context on the codes.

clear

%% 5.2 Converting a PDE to a PIE
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
clear DDE
DDE.A0 = [-1.5, 0; 0.5, -1];
DDE.Adi{1} = [3, 2.25; 0, 0.5]; DDE.tau(1) = 1;
DDE.Adi{2} = [-1, 0; 0, -1];    DDE.tau(2) = 2;
PIE = convert_PIETOOLS_DDE(DDE,'pie')


%% 5.4 Converting an Input-Output system to a PIE
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