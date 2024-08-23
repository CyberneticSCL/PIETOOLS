%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_all_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of computational grid, initial conditions, forcing
% functions and PIE operators in 2D
%
% Inputs:
% 1) PIE - PIE structure of the problem
% 2) uinput - user-defined boundary inputs, forcing and initial conditions
% 3) psize - size of the PIE problem: all variables defining the size of the PIE problem
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% all PIE operators
% 2) coeff - Chebyshev coefficients for initial conditions and forcing functions
% 3) grid - contains physical and computational grid for states differentiable up to order zero (corresponding to a orimary = PDE state discretization)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 04_16_2024

 function [Dop, coeff, grid]=PIESIM_discretize_all_2D(PIE, uinput, psize);

% Discretize computational domain
[grid,gridall]=PIESIM_discretize_domain_2D(uinput,psize);
% n1grid, n2grid are local variables denoting auxiliary computational grid
% of n1, n2 states (not needed outside PIESIM_discretize_all)

% Discretize initial conditions and non-polynomial in space forcing matrix operator

coeff=PIESIM_discretize_icf_2D(uinput,psize,gridall);

% Discretize PIE operators
Dop=PIESIM_discretize_ops_2D(PIE,psize);

 end

