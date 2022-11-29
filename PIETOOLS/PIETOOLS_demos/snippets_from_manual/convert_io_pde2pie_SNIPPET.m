% Code Snippet illustrating usage of convert.m
% See Section 5.2 in Manual for a description.
%
% This document illustrates the process of declaring a PDE with inputs and
% outputs, and converting it to a PIE.

% Equations
% PDE:     \dot(x)(t,s) = 0.5*x_ss(t,s) + s*(2-s)*w(t),         a<s<b
% Outputs: z(t) = int_{a}^{b} x(t,s) ds,
%          y(t) = x(t,1),
% BCs:     x(t,a) = u(t),   
%          x_{s}(t,b) = 0.
% Where w is a disturbance, u a controlled input, z a regulated output, y
% an observed output, and where we use domain:  a=0, b=1.

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
w = state('in');    u = state('in');
z = state('out');   y = state('out');


%%% Second step: Declare the PDE equations and output equations
PDE_dyn = diff(x,t) == 0.5*diff(x,s,2) + s*(2-s)*w;
PDE_z = z == int(x,s,[a,b]);
PDE_y = y == subs(x,s,b);
PDE = addequation(PDE,[PDE_dyn; PDE_z; PDE_y]);
             
                   
%%% Third step: Declare the boundary conditions
PDE_BCs = [subs(x,s,a) == u; 
           subs(diff(x,s),s,b) == 0];
PDE = addequation(PDE,PDE_BCs);


%%% Fourth step: Set the domain of the PDE, and declare observed outputs
%%% and controlled inputs
PDE.dom = [a,b];
PDE = setControl(PDE,u);  PDE = setObserve(PDE,y);

                   
%%% Final step: Convert to a PIE
% This computes the PIE operators T, Tw, Tu, A, B1, B2, C1, D11, D12, C2, 
% D21, D22 defining the PIE representation
PIE = convert(PDE,'pie');


%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off