%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver_PIESIM.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PIESIM Version of PIETOOLS for 1D and 2D problems: 
% For support, contact Y. Peet, Arizona State University at ypeet@asu.edu

% This is the main solver for PIESIM code. It can be called instead
% of PIETOOLS_PDE if a PDE or PDE+ODE soluton is the only required task. 
% Examples can be drawn from 'examples_pde_library_PIESIM.m' for 1D, 
% 'examples_pde_library_PIESIM_2D.m' for 2D, or set up
%  manually within this code (by following the instructions below)


% This routine is the main driver for the solution of the PDE (and PDE+ODE) systems in 1D and 2D in a
% PIE representation framework. The user needs to set up the problem (or call if from examples libary).
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
% PDE/ODE-PDE/DDE

% Set up the dimenion of the problem (dim=1 for 1D problems or dim=2 for 2D problems)

dim=2;

% OPTION 1 (set "init_option=1"): choose one of the examples from the
% library (examples_pde_library_PIESIM) 
% OPTION 2 (set "init_option=2"):
% Define your own input further after init_option=2

init_option=1;
%------------------------------------------------------------------------------
% For 1D problems: example = xxx (between 1 and 38) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_1D.m'
% For 2D problems: example = xxx (between 1 and 6) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_2D.m'
if (init_option==1)
    example=19;

    if (dim==1)
    if (example<1|example>38)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_1D(example);
    else   % dim=2
    if (example<1|example>20)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_2D(example);
    end
end


if (init_option==2)
    syms sx st;
    %------------------------------------------------------------------------------
    % OPTION 2: Place your own input structure below
    %------------------------------------------------------------------------------
    %
    % A required input: PDE, DDE or PIE structure (consult PIETOOLS manual)
    %
    % Optional inputs: simulation options (opts), user-defined inputs
    % (uinput)
    % NOTE: if initial conditions are not specified, they will be
    % defaulted to zero for ODE states and sine waves for distributed
    % states.
%-----------------------------------------------------
    % NOTE: examples below correspond to 1D
%-----------------------------------------------------    
    % DDE example: 
    % tauN=1.3; 
    % DDE.A0=[0 1;-1 .1]; 
    % DDE.Ai{1}=[0 0; -1 0];
    % DDE.Ai{2}=[0 0;1 0]; 
    % DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;
%-----------------------------------------------------
    % PDE example converted to PIE:
%     a=-0.6;b=0.6;
%     PDE.dom=[a b];
%     PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nw=2; PDE.nu=0; PDE.no=0;
%     visc = 0.5;
%     PDE.A0=0; PDE.A1=0; PDE.A2=visc;
%     PDE.B=[1 0 0 0;0 1 0 0];
%     PDE.Bw=[1 0;0 1];
%     
%     uinput.exact =  sin(pi*sx)*exp(-visc*pi^2*st);
%     % Initial conditions for the primary states of the PDE
%     uinput.ic.PDE=  -pi^2*sin(pi*sx);
%     uinput.w(1)= sin(pi*a)*exp(-visc*pi^2*st);
%     uinput.w(2)= sin(pi*b)*exp(-visc*pi^2*st);
%     PIE=convert_PIETOOLS_PDE_batch(PDE);
%     PDE.n.n_pde=[0 0 1];
%--------------------------------------------------------

% Example where left-hand side operator of the PDE dynamics 
% is not the same as T operator (map between fundamental and primary
% states), i.e. PDE equation is not necessarily in the form 
% u_t=Partial derivatives(u), such as below. In this case T operator
% corresponds to the LHS operator, and T0 to the map between the states - 
%
% For more details, see C. Edwards and Y. T. Peet, "Linear Stability of Plane Poiseuille Flow
% in the Sense of Lyapunov", 2023 62nd IEEE Conference on Decision and Control, Marina Bay Sands, Singapore, 4693-4698

% u(x,t)=sin(x^2)*exp(A*t) - solution for u_xxt+alpha*u_t=A*uxx+alpha*A*u

% Boundary conditions: u(a,t)=0, u_x(a,t)=0
            
%    
%              a=0;b=1; 
%              PDE.dom=[a b]; 
%              A=-5;
%              alpha=10;
% 
%              % Batch format
% 
%                PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nw=0; PDE.nu=0; PDE.no=0; 
%                PDE.A0=A; PDE.A1=0; PDE.A2=A/alpha;
%                PDE.B=[1 0 0 0;0 0 1 0];
% 
% %           Exact solution, initial conditions and inhomogeneous inputs   
%              
%               uinput.exact =  sin(sx^2)*exp(A*st);
%              % Initial conditions for the primary states of the PDE
%               uinput.ic.PDE=  diff(sin(sx^2),sx,2);
% 
% 
% 
%         PIE=convert_PIETOOLS_PDE_batch(PDE);
%         PIE.T0=PIE.T;
%         PIE.T.R.R0=PDE.A2/PDE.A0;
%         PDE.n.n_pde=[0,0,1];


end % User-defined initialization option


% If exact solution is known and is desired to be provided for testing,
% select ifexact=true and provide an exact solution as a function in
% symolic form using symbolic variables sx (for space), st (for time) 
% NOTE: Only choose ifexact=true if exact solution for all the states is
%       available, Otherwise, choose ifexact=false.

%exact(1) = sin(pi*sx)*exp(-nu*pi^2*st);
uinput.ifexact=true;
%uinput.ifexact=false;

%-----------------------------------------------------------
opts.plot='yes';
opts.ploteig='yes';

% Here we define parameters related to simulation, do not modify these
% unless necessary 
% Input N - the Chebyshev polynomial discretization order of the
%            distributed states
opts.N=16;

%-----------------------------------------------------------
% Input the desired final time of the solution
opts.tf=0.8;

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
%  schemes, and higher order more accurate schemes Also input the desired
%  time step (dt)

opts.intScheme=1;

if (opts.intScheme==1)
    opts.Norder = 2;
    dt=1e-3;
    % NOTE: if tf is not divisible by dt, dt will be adjusted to a closest
    %       dt value to yield an integer number of time step
    Nsteps=floor(opts.tf/dt);
    dt=opts.tf/Nsteps;
    opts.dt=dt;
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



