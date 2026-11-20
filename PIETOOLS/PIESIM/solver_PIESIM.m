%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver_PIESIM.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PIESIM Version of PIETOOLS for 1D and 2D problems: 
% For support, contact Y. Peet, Arizona State University at ypeet@asu.edu

% This is the wrapper for the PIESIM code if a numerical solution of PDE/ODE problem is the only task required.
% It simulates an example of the 1D or 2D PIESIM example library, chosen by
% the variables dim and example below.
%
% This routine performs:
%
% 1) Conversion of a PDE (or PDE+ODE) problem to a PIE representation
% 2) Spatial discretization of the PIE operators (with Chebyshev polynomial approximation - high-order)
% 3) Temporal discretization (up to 4th order - user-defined parameter) and time advancement of the resulting ODE
% system 
% 4) Plotting and output of results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - solver_PIESIM
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
% If you modify this code, document all changes carefully and include date,
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021
% YP - Added 2D functionality - 4_16_2024

% YP 12/21/2025 - (1) moved options before example setup, (2) pass options to examples for
% generating data array disturbance, (3) moved uinput.ifexact to
% opts.ifexact
% MMP, 11/19/2026: Release comment pass (comment-only).

clear all;
clc
close all;
format long;

%--------------------------------------------------------------
% SETUP OF THE SIMULATIONS USER INPUT BEGINS
%--------------------------------------------------------------

% This script can be used to simulate the following linear systems:
% PDE/ODE or PDE/DDE

% Set up the dimenion of the problem (dim=1 for 1D problems or dim=2 for 2D problems).

dim=1;

if (dim~=1&dim~=2)
    error('Dimension must be specified as 1 or 2');
end

% Set up options for the simulations

% If exact solution is known and is desired to be provided for testing,
% select opts.ifexact=true.
% NOTE: Only choose opts.ifexact=true if exact solution for all the states is
%       available. Otherwise, choose opts.ifexact=false.

opts.ifexact=true;
%-----------------------------------------------------------
opts.plot='yes';
opts.ploteig='yes';

% Here we define parameters related to simulation.

% Input opts.N - the Chebyshev polynomial discretization order of the
%            distributed states 

% For 2D: if different order of discretization is required in different directions,
% enter opts.N as an array, such as opts.N=[16,32]

opts.N=16;

%-----------------------------------------------------------
% Input the desired final time of the solution
opts.tf=0.1;

%-----------------------------------------------------------
% Choose temporal integration scheme 
%  opts.intScheme = 1 - Backward Differentiation Formula (BDF) 
%  opts.intScheme = 2 - Analytical integration in symbolic form 
% Note: opts.intScheme=2 will only work if the boundary and forcing inputs
%       are simple integrable functions of time, and the matrix 
%       Atotal=inv(T)*A is diagonalizable. An error will be issued if 
%       matrix is not diagonalizable, and a default integration scheme 
%       given by opts.intScheme = 1 (BDF) of order 2 (opts.Norder=2) will 
%       be executed.

% Choose opts.intScheme=1 if a
%  temporal history of solution is required (solution history is not
%  provided with analytical integration)
% if opts.intScheme = 1
%  Choose the order of numerical time integration scheme (Norder). A time
%  integration scheme available is Backward Differentiation Formula (BDF).
%  opts.Norder = 1, 2, 3 or 4 can be chosen. Lower order yield more robust
%  schemes, and higher order more accurate schemes. Also input the desired
%  time step (opts.dt)

opts.intScheme=1;

if (opts.intScheme==1)
    opts.Norder = 2;
    opts.dt=1e-3;
end


% Set up example to run from the examples library (additional examples can
% be added to the library by the user).
%------------------------------------------------------------------------------
% For 1D problems: example = xxx (between 1 and 41) to correspond to an Example number in
% the 1D example library
% For 2D problems: example = xxx (between 1 and 32) to correspond to an Example number in
% the 2D example library
%------------------------------------------------------------------------------

    example=41;

    if (dim==1)
    if (example<1|example>41)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_1D(example);
    else   % dim==2
    if (example<1|example>32)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_2D(example);
    end

%--------------------------------------------------------------------------
%   USER INPUT ENDS
%--------------------------------------------------------------------------

% If disturbances are used while supplying PIE input to PIESIM, 
% user needs to ensure that finite-dimensional distubrances are
% entered before infinite-dimensional disturbances for a proper
% functionality. Disturbances are not reordered automatically with the PIE
% structure as they are with the PDE structure.

if exist('PIE','var')
    solution = PIESIM(PIE,opts,uinput);
elseif exist('DDE','var')
    solution=PIESIM(DDE,opts,uinput);
else
    [solution,grid] = PIESIM(PDE,opts,uinput);
end

