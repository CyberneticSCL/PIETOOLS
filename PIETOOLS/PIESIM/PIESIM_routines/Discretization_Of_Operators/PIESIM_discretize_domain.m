%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_domain.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of the computational domain  
% gridall{i} holds the N-i+2 Chebyshev-Gauss-Lobatto nodes
% xi_k = cos(pi*k/(N-i+1)), k = 0,...,N-i+1, mapped to [a,b] by
% s = (b-a)/2*xi + (b+a)/2: the grid of the states differentiable i-1 times.
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
%    ('struct' with fields a, b, the domain [a,b])
% 2) psize - size of the problem
%    ('struct' with fields N, the Chebyshev order, and n, n(i) the number
%    of states differentiable i-1 times)
%
% Outputs:
% gridall -cell array containing grids of different degrees of differentiability
%    (1 x length(psize.n) cell; gridall{i} (N-i+2) x 1 'double')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_domain
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
% Initial coding YP  - 11_05_2021
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output grid removed;
%   help and comments revised.
function gridall=PIESIM_discretize_domain(uinput,psize);

% Define local variables

a=uinput.a;
b=uinput.b;
N=psize.N;


% Define grids in the physical space as vector arrays. Spatial grid points
% are colocated with the Chebyshev nodes.

% Define Chebyshev grid in computational domain for each type of state
% Then convert it to physical domain

for i=1:length(psize.n)
grid_comp = cos(pi*(0:N-i+1)/(N-i+1))';
gridall{i}=0.5*(b-a)*grid_comp+0.5*(b+a);
end