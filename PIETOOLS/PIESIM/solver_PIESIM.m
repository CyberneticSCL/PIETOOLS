%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver_PIESIM.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PIESIM Version of PIETOOLS for 1D and 2D problems: 
% For support, contact Y. Peet, Arizona State University at ypeet@asu.edu

% This is the main driver for the PIESIM code if a numerical solution of PDE/ODE problem is the only task required.
% It can be called instead of PIETOOLS_PDE in this case. 
% Examples can be drawn from 'examples_pde_library_PIESIM.m' for 1D or
% 'examples_pde_library_PIESIM_2D.m' for 2D.
%
% This routine performs:
%
% 1) Conversion of a PDE (or PDE+ODE) problem to a PIE representation
% 2) Spatial discretization of the PIE operators (with Chebyshev polynomial approximation - high-order)
% 3) Temporal discretization (up to 4th order - user-defined parameter) and time advancement of the resulting ODE
% system 
% 4) Plotting and output of results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date,
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021
% YP - Added 2D functionality - 4_16_2024

clear all;
clc;
close all;
format long;

% Set up the dimension of the problem 

%--------------------------------------------------------------
% SETUP OF THE SIMULATIONS USER INPUT BEGINS
%--------------------------------------------------------------

% This script can be used to simulate the following linear systems:
% PDE/ODE or PDE/DDE

% Set up the dimenion of the problem (dim=1 for 1D problems or dim=2 for 2D problems).

dim=1;

% Set up example to run from the examples library (additional examples can
% be added to the library by the user).
%------------------------------------------------------------------------------
% For 1D problems: example = xxx (between 1 and 38) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_1D.m'
% For 2D problems: example = xxx (between 1 and 19) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_2D.m'
%------------------------------------------------------------------------------

    example=1;

    if (dim==1)
    if (example<1|example>38)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_1D(example);
    else   % dim=2
    if (example<1|example>19)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_2D(example);
    end


% If exact solution is known (see examles) and is desired to be provided for testing,
% select uinput.ifexact=true.
% NOTE: Only choose uinput.ifexact=true if exact solution for all the states is
%       available, Otherwise, choose ifexact=false.

uinput.ifexact=true;
%-----------------------------------------------------------
opts.plot='yes';
opts.ploteig='yes';

% Here we define parameters related to simulation.

% Input N - the Chebyshev polynomial discretization order of the
%            distributed states
opts.N=16;

%-----------------------------------------------------------
% Input the desired final time of the solution
opts.tf=1.6;

%-----------------------------------------------------------
% Choose temporal integration scheme 
%  opts.intScheme = 1 - Backward Differentiation Formula (BDF) 
%  opts.intScheme = 2 - Analytical integration in symbolic form 
% Note: opts.intScheme=2 will only work if the boundary and forcing inputs
%       are simple integrable functions of time, and the matrix 
%       Atotal=inv(M)*A is diagonalizable. An error will be issued if 
%       matrix is not diagonalizable, and a default integration scheme 
%       given by opts.intScheme = 1 (BDF) of order 2 (opts.Norder=2) will 
%       be executed.

% Choose opts.intScheme=1 if a
%  temporal history of solution is required (solution history is not
%  provided with analytical integration)
% if opts.intScheme = 1
%  Choose the order of numerical time integration scheme (Norder). A time
%  integration scheme available is Backward Differentiation Formula (BDF).
%  Norder = 1, 2, 3 or 4 can be chosen. Lower order yield more robust
%  schemes, and higher order more accurate schemes. Also input the desired
%  time step (dt)

opts.intScheme=1;

if (opts.intScheme==1)
    opts.Norder = 2;
    opts.dt=0.02;
end

%--------------------------------------------------------------------------
%   USER INPUT ENDS
%--------------------------------------------------------------------------

if exist('PIE','var')
    solution = PIESIM(PIE,opts,uinput,PDE.n.n_pde);
elseif exist('DDE','var')
    solution=PIESIM(DDE,opts,uinput);
else
    solution = PIESIM(PDE,opts,uinput);
end



