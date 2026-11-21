%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_all_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of computational grid, initial conditions, forcing
% functions and PIE operators in 2D
%
% Inputs:
% 1) PIE - PIE structure of the problem
%    ('struct' with the operators T, Tw, Tu, A, B1, B2, C1, C2, D11, D12,
%    D21, D22, each a 'struct' of its kernels R00, R0x, ..., R22 with
%    fields dim, var1, var2, on [-1,1]^2; and misc)
% 2) uinput - user-defined boundary inputs, forcing and initial conditions
%    ('struct' with fields dom (2x2, the domain [a,b]x[c,d]), ic, w, u,
%    wspace, uspace, wsep, usep)
% 3) psize - size of the PIE problem: all variables defining the size of the PIE problem
%    ('struct'; N = [N_x, N_y] the Chebyshev orders, n0, nx, ny, n the
%    numbers of ODE states and of states on [a,b], [c,d] and
%    [a,b]x[c,d] by order of differentiability, nw, nu (nw0, nwx, nwy, nw2,
%    nu0, nux, nuy, nu2) the numbers of disturbances and control inputs)
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% all PIE operators
%    ('struct')
% 2) coeff - Chebyshev coefficients for initial conditions and forcing functions
%    ('struct' with fields acheb_f0, w, u)
% 3) gridall - physical grids for solution states of different degree of
% differentiability
%    ('struct' with cell fields x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_all_2D
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
% Initial coding YP  - 04_16_2024
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: third output grid
%   replaced by gridall, the physical grids of the states of every degree of
%   differentiability; help and comments revised.

 function [Dop, coeff, gridall]=PIESIM_discretize_all_2D(PIE, uinput, psize);

% Discretize computational domain
gridall=PIESIM_discretize_domain_2D(uinput,psize);

% Discretize initial conditions and non-polynomial in space forcing matrix operator

coeff=PIESIM_discretize_icf_2D(uinput,psize,gridall);

% Discretize PIE operators
Dop=PIESIM_discretize_ops_2D(PIE,psize);

 end

