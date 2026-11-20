%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_domain_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of the computational domain in 2D  
% gridall.x{i} holds the N_x-i+2 Chebyshev-Gauss-Lobatto nodes
% cos(pi*k/(N_x-i+1)), k = 0,...,N_x-i+1, mapped to [a,b] by
% x = (b-a)/2*xi + (b+a)/2; gridall.y{i} likewise with N_y on [c,d]: the
% grids of the states differentiable i-1 times in that direction.
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
%    ('struct' with field dom, 2x2 'double' [a b; c d])
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
%    ('struct' with fields N = [N_x, N_y], the Chebyshev orders, and n, nx,
%    ny, the numbers of states on [a,b]x[c,d], [a,b], [c,d] by order of
%    differentiability)
%
% Outputs:
% 1) gridall - field containing physical grids of different degree of
% differentiabiloity
%  gridall.x - cell array containing grids in x direction of different degrees of differentiability
%  gridall.y - cell array containing grids in y direction of different degrees of differentiability
%  (i = 1,...,max(length(n), length(nx), length(ny)); x{i} (N_x-i+2) x 1,
%  y{i} (N_y-i+2) x 1 'double')
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_domain_2D
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
% Initial coding YP  - 04/16/2024
% 2/12/2026: YP - modified grids handling to account for different degree
% of differentiability and number of grid points
% MMP, 11/19/2026: Release comment pass (comment-only).

function gridall=PIESIM_discretize_domain_2D(uinput,psize);

% Define local variables

a=uinput.dom(1,1);
b=uinput.dom(1,2);
c=uinput.dom(2,1);
d=uinput.dom(2,2);
N=psize.N;

% Define grids in the physical space as vector arrays. Spatial grid points
% are colocated with the Chebyshev nodes.

% Define Chebyshev grid in computational domain for each type of state
% Then convert it to physical domain


for i=1:max([length(psize.n),length(psize.nx),length(psize.ny)])
grid_comp_x = cos(pi*(0:N(1)-i+1)/(N(1)-i+1))';
grid_comp_y = cos(pi*(0:N(2)-i+1)/(N(2)-i+1))';
gridall.x{i}=0.5*(b-a)*grid_comp_x+0.5*(b+a);
gridall.y{i}=0.5*(d-c)*grid_comp_y+0.5*(d+c);
end

