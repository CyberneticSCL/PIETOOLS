%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_all.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of computational grid, initial conditions, forcing functions and PIE operators 
%
% Inputs:
% 1) PIE - PIE structure of the problem
% 2) uinput - user-defined boundary inputs, forcing and initial conditions
% 3) psize - size of the problem: contains the variables nu,nw,no,N,n0,n1,n2
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% T,Tu,Tw,A,Bi,Ci
% 2) coeff - Chebyshev coefficients for initial conditions and forcing functions
% 3) grid - contains physical and computational grid for n0 states
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_18_2021
function [Dop, coeff, grid]=PIESIM_discretize_all(PIE, uinput, psize);

% Discretize computational domain
[grid,gridall]=PIESIM_discretize_domain(uinput,psize);
% n1grid, n2grid are local variables denoting auxiliary computational grid
% of n1, n2 states (not needed outside PIESIM_discretize_all)

% Discretize initial conditions and non-polynomial in space forcing matrix operator

[coeff, B1_nonpol]=PIESIM_discretize_icf(uinput,psize,gridall);

% Discretize PIE operators
Dop=PIESIM_discretize_ops(PIE,psize);

% Add B1 contrbution from non-polynomial in space forcing


if ~isempty(B1_nonpol) 
    Dop.B1cheb=Dop.B1cheb+B1_nonpol;
end
