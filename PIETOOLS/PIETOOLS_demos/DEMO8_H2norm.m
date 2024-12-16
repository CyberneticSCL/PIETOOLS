%%
clear all; clc;close all;
%%
% DEMO_Simple_H2Norm_Problem.m
% See the manual for full description.
%
% This document illustrates, with a simple example, how PIETOOLS can be used to
% analyse and design optimal controllers and observers based on the H2 norm
%  The example is a pure convection equation in 1D: 
%  PDE :         x_{t} = x_{s} + (s-s^2)w(t)
%  With BC     x(s=1) = 0
%  And output  z(t) = int(x(t,s),s,0,1)
%% First, declare the independent variables as polynomial variables
pvar t s;
%, creat the structures for state, input and output,
x=state('pde');
w=state('in');
z=state('out');
% initialize the pde structure,
pde = sys();
% declare the dynamic equation,
eq_dyn=diff(x,t,1)==diff(x,s,1)+(s-s^2)*w;
% declare the output equation,
eq_out=z==int(x,s,[0,1]);
% declare the boundary condition
bc = subs(x,s,1)==0;
% add all system equations to the pde structure,
pde = addequation(pde,[eq_dyn;eq_out;bc]);
% and finally converts the pde to PIE representation
PIE = convert(pde,'pie');   
C1=PIE.C1; 
A=PIE.A;
B1=PIE.B1;
T=PIE.T;
% Initialize the program structure
prog = lpiprogram(PIE.vars,PIE.dom); 
% declare the polynomial decision variable,
dpvar gam;
% sets gam as decision variable of the sosprogram,
prog = lpidecvar(prog, gam);
% Creats positive PI decision variable W : W(s) \in L_2 with
% default polynomial degrees up to 3.
[prog, W] = poslpivar(prog, [0,1]);
% Declare the LPI
Dop =  A*W*T'+T*W*A'+B1*B1';
% Impose positiviness on -Dop [1], adding the first inequality to the sosprogram.
% By default the operators are choosen positive on any domain [a,b].
prog = lpi_ineq(prog,-Dop);
Aux=C1*W*C1';
traceVal= trace(Aux.P);
% Adds the second inequality to the sosprogram
prog = lpi_ineq(prog, gam-traceVal);
% Finally, sets the solver,
sos_opts.solver='mosek';
%sets gam as the objective function of the minimization problem
prog = lpisetobj(prog, gam); 
% and search for a polynomial solution
prog = lpisolve(prog,sos_opts); 
% If primal and dual feasible, retrieve the controlability grammian and the
% upper bound on the H2 norm
Wc = getsol_lpivar(prog,W);
gamd = sqrt(double(lpigetsol(prog,gam)));
%% [1] M. Peet, “A partial integral equation (PIE) representation of
% coupled linear PDEs and scalable stability analysis using LMIs,”
% Automatica, vol. 125, p. 109473, 2021.