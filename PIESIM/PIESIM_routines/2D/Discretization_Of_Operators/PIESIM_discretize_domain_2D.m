%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_domain_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of the computational domain in 2D  
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
%
% Outputs:
% 1) grid - physical and computational grid for states differentiable up to order zero (corresponding to a primary = PDE state discretization)
% 2) gridall - cell array of size dmax containing physical grid for all states
% depending on their degree of differentiability; dmax corresponds to the
% maximum degree of differentiability among all the states
%  gridall.x - grids in x direction
%  gridall.y - grids in y direction
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 16_04_2024

function [grid,gridall]=PIESIM_discretize_domain_2D(uinput,psize);

% Define local variables

a=uinput.dom(1,1);
b=uinput.dom(1,2);
c=uinput.dom(2,1);
d=uinput.dom(2,2);
N=psize.N;

% Assume for now that state differentiability is the same in x and y. Seems
% like what it's done in Declan't paper
% Define grids in the physical space as vector arrays. Spatial grid points
% are colocated with the Chebyshev nodes.

% Define Chebyshev grid in computational domain for each type of state
% Then convert it to physical domain


for i=1:max([length(psize.n),length(psize.nx),length(psize.ny)])
grid_comp{i} = cos(pi*(0:N-i+1)/(N-i+1))';
gridall.x{i}=0.5*(b-a)*grid_comp{i}+0.5*(b+a);
gridall.y{i}=0.5*(d-c)*grid_comp{i}+0.5*(d+c);
end

grid.phys=[gridall.x{1},gridall.y{1}];
grid.comp=grid_comp{1};