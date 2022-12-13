% Code Snippet illustrating usage of convert.m
% See Section 5.2 in Manual for a description.
%
% This document illustrates the process of declaring a heat equation
% with integral BCs, and converting it to a PIE.

% Equations
% PDE:     \dot(x)(t,s) = x_ss(t,s),         a<s<b
% BCs:     x(t,a) + int(x,s,[a,b]) = 0,   
%          x(t,b) + int(x,s,[a,b]) = 0.
% Domain:  a=0, b=1

%%
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%% Zeroth step: Define temporal variable t and spatial variable s in [a,b]
pvar s t
a = 0;   b = 1;


%%% First step: Declare the state components, inputs and outputs (if
%%% applicable)
PDE = sys();
x = state('pde');


%%% Second step: Declare the PDE equations
PDE_dyn = diff(x,t) == diff(x,s,2);
PDE = addequation(PDE,PDE_dyn);
             
                   
%%% Third step: Declare the boundary conditions
PDE_BCs = [subs(x,s,0) + int(x,s,[a,b]) == 0;
           subs(x,s,1) + int(x,s,[a,b]) == 0];
PDE = addequation(PDE,PDE_BCs);


%%% Fourth step: Set the domain of the PDE
PDE.dom = [a,b];
         
                   
%%% Final step: Convert to a PIE
% This computes the PIE operators T, Tw, Tu, A, B1, B2, C1, D11, D12, C2, 
% D21, D22 defining the PIE representation
PIE = convert(PDE,'pie');


%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off