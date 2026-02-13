%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_domain_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of the computational domain in 2D  
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
%
% Outputs:
% 1) grid - field containing the following sub-fields:
% grid.phys - physical grid for states differentiable up to order zero (corresponding to a primary = PDE state discretization)
%  grid.x - cell array containing grids in x direction of different degrees of differentiability
%  grid.y - cell array containing grids in y direction of different degrees of differentiability
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 16_04_2024
% 2/12/2026: YP - modified grids handling

function grid=PIESIM_discretize_domain_2D(uinput,psize);

% Define local variables

a=uinput.dom(1,1);
b=uinput.dom(1,2);
c=uinput.dom(2,1);
d=uinput.dom(2,2);
N=psize.N;

% Assume for now that state differentiability is the same in x and y. Seems
% like what it's done in Declan's paper
% Define grids in the physical space as vector arrays. Spatial grid points
% are colocated with the Chebyshev nodes.

% Define Chebyshev grid in computational domain for each type of state
% Then convert it to physical domain


for i=1:max([length(psize.n),length(psize.nx),length(psize.ny)])
grid_comp_x = cos(pi*(0:N(1)-i+1)/(N(1)-i+1))';
grid_comp_y = cos(pi*(0:N(2)-i+1)/(N(2)-i+1))';
<<<<<<< HEAD
grid.x{i}=0.5*(b-a)*grid_comp_x+0.5*(b+a);
grid.y{i}=0.5*(d-c)*grid_comp_y+0.5*(d+c);
=======
gridall.x{i}=0.5*(b-a)*grid_comp_x+0.5*(b+a);
gridall.y{i}=0.5*(d-c)*grid_comp_y+0.5*(d+c);
>>>>>>> e5f7b8c94688f4b02e0f00f92dfa09490ed2bb61
end

grid.phys={grid.x{1};grid.y{1}};

