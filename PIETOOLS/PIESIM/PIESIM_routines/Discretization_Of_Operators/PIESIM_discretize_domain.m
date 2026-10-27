%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_domain.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of the computational domain  
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
% 2) psize - size of the problem
%
% Outputs:
% gridall -cell array containing grids of different degrees of differentiability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 11_05_2021
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