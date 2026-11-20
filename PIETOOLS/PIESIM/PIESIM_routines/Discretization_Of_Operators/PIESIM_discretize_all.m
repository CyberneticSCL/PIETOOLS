%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_all.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of computational grid, initial conditions, forcing functions and PIE operators 
%
% Inputs:
% 1) PIE - PIE structure of the problem
%    ('struct' with the operators T, Tw, Tu, A, B1, B2, C1, C2, D11, D12,
%    D21, D22, each a 'struct' of its kernels P, Q1, Q2, R.R0, R.R1, R.R2
%    with the field dim, on [-1,1]; and misc)
% 2) uinput - user-defined boundary inputs, forcing and initial conditions
%    ('struct' with fields a, b (the domain [a,b]), ic, wspace, uspace and,
%    optionally, Bpw_nonpol (non-polynomial forcing))
% 3) psize - size of the problem: contains the variables nu,nw,N,n0,n
%    ('struct'; N the Chebyshev order, n0 the number of ODE states, n(i)
%    the number of PDE states differentiable i-1 times, nw, nu (nw0, nwx,
%    nu0, nux) the numbers of disturbances and control inputs (finite-
%    dimensional, functions of s))
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% T,Tu,Tw,A,Bi,Ci
%    ('struct'; B1cheb includes the non-polynomial forcing)
% 2) coeff - Chebyshev coefficients for initial conditions and forcing functions
%    ('struct' with fields acheb_f0, w, u)
% 4) gridall - contains all grids for all states with different degree of
% differentiability (internal)
%    (1 x length(psize.n) cell; gridall{i} the grid of the states
%    differentiable i-1 times)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_all
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 9_18_2021
% MMP, 11/19/2026: Release comment pass (comment-only).
function [Dop, coeff, gridall]=PIESIM_discretize_all(PIE, uinput, psize);

% Discretize computational domain
gridall=PIESIM_discretize_domain(uinput,psize);

% Discretize initial conditions and non-polynomial in space forcing matrix operator

[coeff, B1_nonpol]=PIESIM_discretize_icf(uinput,psize,gridall);

% Discretize PIE operators
Dop=PIESIM_discretize_ops(PIE,psize);

% Add B1 contrbution from non-polynomial in space forcing


if ~isempty(B1_nonpol) 
    Dop.B1cheb=Dop.B1cheb+B1_nonpol;
end
