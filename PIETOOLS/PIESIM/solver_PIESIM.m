%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver_PIESIM.m     PIETOOLS 2021d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PIESIM Version: 2021d
% date: 9/6/2021 
% NOTE:
% For support, contact Y. Peet, Arizona State University at ypeet@asu.edu


% This is the main solver for PIESIM code. 
% It perfoms the following functions:
% 1) Sets up the PDE problem (either by calling an example from the library
% or setting it up manually)
% 2) Call the executive solver for PIESIM (exectuive_PIESIM.m)

% The executive rotuine performs the following functions:
% 1) Checks if all inputs are defined and (if not) setting up default
% options
% 1) Converts the PDE to PIE problem
% 2) Discretizes the PIE operators with Chebyshev polynomials
% 3) Temporally integrates the corresponding spatially-discretized
% equations (3 temporal schemes are possible: BDF, Gauss integration, and
% analytical integraiton - see options in solver_PIESIM.m)
% 4) Outputs and plots the PDE and ODE solutions at a final time and
% (if opts.intScheme=1) - plots the ODE solutions versus time
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date,
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021
% 
 clear;
 clc;
 close all;
format long;

% addpath(genpath('.')) % makes sure any local version  of SOSTOOLS in the current folder is at the head of the path

% System definition: a system of second-order PDEs 

% \dot [x(t,s)] =  A0(s) [x_0(t,s)] + A1(s) [x_1s(t,s)] + A2(s) [x_2ss(t,s)]+ f(t,s)

% NOTE: only PDEs with polynomial coefficients are currently supported  -  that is, matrices A0, A1
% and A2 only admit polynomial entries 

%--------------------------------------------------------------
% SETUP OF THE SIMULATIONS
% USER INPUT BEGINS
%--------------------------------------------------------------

% Two options to do initialization of the PDE/PDE+ODE/DDE problem.
 
 % OPTION 1 (set "init_option=1"): choose one of the examples from the library (examples_pde_library_PIESIM)
 % example = xxx (between 1 and 33) to correspond to an Example number in the Library 
 
 
 % OPTION 2 (set "init_option=2"): Define your own input structure below 
 
 init_option=1;
 %------------------------------------------------------------------------------
 
 % OPTION 1: choose one of the examples from the library (examples_pde_library_PIESIM)
 % example = xxx (between 1 and 33) to correspond to an Example number in the Library 
 
 if (init_option==1)
 example=2;
 
 if (example<1|example>33)
 disp('Warning: Example number is outside of the range. Defaulting to example=1');
 example=1;
 end
 [PDE,uinput]=examples_pde_library_PIESIM(example);
 end

 
 if (init_option==2)
%------------------------------------------------------------------------------
% OPTION 2: Place your own input structure below
%------------------------------------------------------------------------------
%
% A required input: PDE structure (consult PIETOOLS manual)
% 
% Optional inputs: initial conditions, inhomogeneous inputs uinput, winput
% NOTE: if initial conditions are not specified, they will be defaulted to
% zero.
%------------------------------------------------------------------------------
%   PUT YOUR INPUT BELOW
%------------------------------------------------------------------------------
%   USER INPUT BEGINS
%------------------------------------------------------------------------------


% tauN=1.3;
% DDE.A0=[0 1;-1 .1];
% DDE.Ai{1}=[0 0; -1 0];
% DDE.Ai{2}=[0 0;1 0];
% DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;

 end % User-defined initialization option

          
 % If exact solution is known and is desired to be provided for testing,
 % select ifexact=true and provide an exact solution as a function in
 % symolic form using symbolic variables sx (for space), st (for time)
 % NOTE: only choose ifexact=true if exact solution for all the states is
 % available. Otherwise, choose ifexact=false.

 %exact(1) = sin(pi*sx)*exp(-nu*pi^2*st);
 
 uinput.ifexact=true;
 
 %uinput.ifexact=false;
 
 opts.plot='yes';
 
 
%------------------------------------------------------------------------------
%   USER INPUT ENDS
%------------------------------------------------------------------------------
 %------------------------------------------------------------

 % Input N - the Chebyshev polynomial discretization order of the primary
 % state
 
 %-------------------------------------------------------------
 
 opts.N=16; 
 
 %-----------------------------------------------------------
 
   % Input the desired final time of the solution
  
  opts.tf=0.1;
 
%  Choose temporal integration scheme
%  opts.intScheme = 1 - Backward Differentiation Formula (BDF)
%  opts.intScheme = 2 - Gauss evaluation of the analytical solution
%  opts.intScheme = 3 - Analytical integration in symbolic form
%  Note: opts.intScheme=3 will only work if the boundary and forcing inputs are simple integrable
%  functions of time, and the matrix Atotal=inv(M)*A is diagonalizable. 
%  An error will be issued if matrix is not diagonalizable, and a default
%  integration scheme given by  opts.intScheme = 1
%  (BDF) of order 2 (opts.Norder=2) will be executed.
%  Note: symbolic integration is also the slowest (however, the most accurate when it works).
%  BDF is the fastest, but less accurate 
%  Choose opts.intScheme=1 if a temporal history of solution is required
%  (this is the only scheme that provides it)

%------------------------------------------------
   opts.intScheme=1;
%------------------------------------------------   
   
 
 % if opts.intScheme = 1
 % Choose the order of numerical time integration scheme (Norder). 
 % A time integration scheme available is Backward Differentiation Formula (BDF).
 % Norder = 1, 2, 3 or 4 can be chosen. Lower order yield more robust schemes, and higher
 % order more accurate schemes
 % Also input the desired time step (dt)

  if (opts.intScheme==1)
  Norder = 2;
  dt=1e-3;
    
  % NOTE: if tf is not divisible by dt, dt will be adjusted to a closest dt
  % value to yield an integer number of time steps
  
  Nsteps=floor(opts.tf/dt);
  dt=opts.tf/Nsteps;
  
  opts.Norder=Norder;
  opts.dt=dt;
  opts.Nsteps=Nsteps;
  end
  
  % if opts.intScheme = 2
  % Choose the number of untervals for Gauss evaluation of the integral (Nint)
  % Choose the number of Gauss integration points per interval (Npgauss)
  
  if (opts.intScheme==2)
      opts.Nint=1;
      opts.Npgauss=100;
  end
  
  
% Check if PDE, opts and input structures are defined
 %if ~exist('PDE','var')
 %error('PDE structure is undefined. It is necessary for executing the code');
 %end
 if ~exist('uinput','var')
 disp('Warning: user input structure ``uinput'' is not defined. Defaulting to zero');
 uinput=struct;
 end

  
 if exist('PDE','var') 
     type=PDE;
 elseif exist('DDE','var') 
     type=DDE;
 end
 
  solution=executive_PIESIM(type,opts,uinput);