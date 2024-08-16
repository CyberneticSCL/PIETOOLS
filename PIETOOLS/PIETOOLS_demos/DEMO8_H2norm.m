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
%%
pvar t s;
x=state('pde');
w=state('in');
z=state('out');
pde = sys();
eq_dyn=diff(x,t,1)==diff(x,s,1)+(s-s^2)*w;
eq_out=z==int(x,s,[0,1]);
pde = addequation(pde,[eq_dyn;eq_out]);
bc = subs(x,s,1)==0;
pde = addequation(pde,bc);
PIE = convert(pde,'pie');   
varlist = [PIE.A.var1; PIE.A.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist); % Initialize the program structure
dpvar gam;
prog = sosdecvar(prog, gam); %this sets gam as decision var
prog = sosineq(prog, gam); %this ensures gamma is lower bounded
prog = sossetobj(prog, gam); %this minimizes gamma
options1.sep = 1; %this is to select separable case, R1=R2
options1.exclude=[0 0 0 0]; %don't exclude any part of the PI operator
[prog, Wop] = poslpivar(prog, [PIE.A.dim(1,1) ,PIE.A.dim(2,1)],PIE.A.I,{[5],[0 0 0],[0 0 0]},options1);
Dop =  (PIE.A*Wop)*PIE.T'+PIE.T*(Wop*PIE.A')+PIE.B1*PIE.B1';
opts.psatz = 0;opts.pure = 1;
prog = lpi_ineq(prog,-Dop,opts);
tempObj = PIE.C1*Wop*PIE.C1';
tempMat = tempObj.P;
traceVal=0;
for idx = 1:size(tempMat,1)
    traceVal = traceVal+tempMat(idx,idx);
end
prog = sosineq(prog, gam-traceVal);
sos_opts.solver='mosek';
prog = sossolve(prog,sos_opts); 
Wc = getsol_lpivar(prog,Wop);
gamd = sqrt(double(sosgetsol(prog,gam)));