% This document has an illustration of converting a reaction-diffusion pde
% equation coupled with an ODE at the boundary to a PIE

% Equations
% ODE:     \dot(xo)(t) = -2xo(t)
% PDE:     \dot(xp)(t,s) = lamb xp(t,s) + xp_ss(t,s),   a<s<b
% BCs:     xp(t,a) =0,  xp(t,b) = xo(t)
% Domain:  a=0, b=1

%%
%%%%%%%%%%%%%%%%%% Code Snippet %%%%%%%%%%%%%%%%%%
clc; clear;

%%% Zero step: Define the polynomial variable s, if needed
pvar s theta;

%%% First step: Define the Domain [a,b]
PDE.dom = [0,1];

%%% Second step: xp is twice differentiated and hence is in n2
% number of pde states are defined as n2=1 and n0=n1=0
PDE.n0=0; PDE.n1=0; PDE.n2=1;

%%% Define PDE parameters
lamb = 5;
PDE.A0 = lamb; PDE.A2 = 1;

%%% Define Boundary conditions in the form B[xp(a); xp(b); xp_s(a) xp_s(b)]= Bxo xo
PDE.B = [1 0 0 0; 0 1 0 0]; 

%%% Use conversion script to convert 
% This defines the PIE operators Top, TB1op, TB2op, Aop, B1op, B2op, C1op,
% D11op, D12op, C2op, D21op, D22op
PIE = convert_PIETOOLS_PDE_batch(PDE);


