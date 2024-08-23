% Description of how to run PIESIM 
%
% NOTES:
% For support, contact Y. Peet, Arizona State University at ypeet@asu.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM - README
%
% Copyright (C)2024  Y. Peet, M. Peet

%----------------------------------------------------------
%
% PIESIM is a simulator based on PIETOOLS - it currently supports 1D and 2D problems. 
% It uses Chebyshev polynomial dicretization in space and several integration options in time

% To use: run `solver_PIESIM.m'
% NOTE: PIESIM can also be run out of PIETOOLS_PDE after solving
% stability, estimation or a controller synthesis problem

% Specify the dimension of the problem: 
% dim=1 for 1D problems or dim=2 for 2D problems

% An extensive library of examples is provided for 1D and 2D problems in 
% the files `examples_pde_library_PIESIM_1D.m' (for
% 1D) or examples_pde_library_PIESIM_2D.m' (for 2D)

% If user desires to provie their own PDE problem, this can be done from
% the PIETOOLS GUI or a command line infrastructure (see User's manual). In
% this case, PIESIM is called out of PIETOOLS_PDE.

% User can also add their own examples to the exampls library. In this
% case, append your example at the end of the file between the last example and the keyword
% 'end'.

% Start your example with the line
% case #, e.g. case 40 (make sure this case number does not exist yet)

% Input that number in the `solver_PIESIM.m' under example = (such as, example = 40).

