% Code Snippet illustrating usage of convert_PIETOOLS_DDE.m
% See Section 5.3 in Manual for a description.
%
% This document illustrates the process of declaring a DDE, 
% and converting it to a PIE.

% Equations
% DDE:  \dot{x}(t) = [-1.5,  0] x(t) + int_{-1}^{0} [3, 2.25] x(t+s) ds + int_{-2}^{0} [-1,  0] x(t+s) ds    
%                    [ 0.5, -1]                     [0, 0.5 ]                          [ 0, -1]  


%%
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%% First step: Declare the DDE
DDE.A0 = [-1.5, 0; 0.5, -1];                        % Non-delayed term  A0*x(t)
DDE.Adi{1} = [3, 2.25; 0, 0.5];     DDE.tau(1) = 1; % Distributed delay int_{-1}^{0} Ad1*x(t+s) ds
DDE.Adi{2} = [-1, 0; 0, -1];        DDE.tau(2) = 2; % Distributed delay int_{-2}^{0} Ad1*x(t+s) ds


%%% Final step: Convert to a PIE
% This computes the PIE operators T, Tw, Tu, A, B1, B2, C1, D11, D12, C2, 
% D21, D22 defining the PIE representation
PIE = convert_PIETOOLS_DDE(DDE,'pie');


%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off