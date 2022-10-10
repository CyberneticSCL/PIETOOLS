% This document illustrates the process of declaring a reaction-diffusion 
% PDE coupled with an ODE at the boundary, and converting it to a PIE.

% Equations
% ODE:     \dot(xo)(t) = -2xo(t)
% PDE:     \dot(xp)(t,s) = lamb xp(t,s) + xp_ss(t,s),   a<s<b
% BCs:     xp(t,a) = 0,   xp(t,b) = xo(t)
% Domain:  a=0, b=1

%%
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%% Zeroth step: Define temporal variable t and spatial variable s in [a,b]
lamb = 5;
pvar t s theta;
a = 0;      b = 1;


%%% First step: Declare the state components, inputs and outputs (if
%%% applicable)
PDE = sys();
xo = state('ode');
xp = state('pde');


%%% Second step: Declare the PDE equations
PDE = addequation(PDE,[diff(xo,t) + 2*xo;
                       diff(xp,t) - lamb*xp - diff(xp,s,2)]);
             
                   
%%% Third step: Declare the boundary conditions
PDE = addequation(PDE,[subs(xp,s,a);
                       subs(xp,s,b) - xo]);
         
                   
%%% Final step: Convert to a PIE
% This computes the PIE operators T, Tw, Tu, A, B1, B2, C1, D11, D12, C2, 
% D21, D22 defining the PIE representation
PIE = convert(PDE,'pie');


%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off