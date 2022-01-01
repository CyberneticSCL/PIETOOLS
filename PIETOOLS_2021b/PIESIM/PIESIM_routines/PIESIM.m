%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please, contact Y. Peet at ypeet@asu.edu for support

% This matlab function performs the following operations:
% 1) Checks if all inputs are defined and (if not) setting up default
% options
% 2) If the input structure is PDE or DDE: converts it to PIE problem
% 3) Discretizes the PIE operators with Chebyshev polynomials
% 4) Temporally integrates the corresponding spatially-discretized
%    equations (3 temporal schemes are possible: BDF, Gauss integration, 
%    and analytical integraiton - see options in solver_PIESIM.m)
% 5) If opts.plot='yes', outputs and plots the PDE/DDE and ODE solutions at a final time and
%    (if opts.intScheme=1) - plots the ODE solutions versus time
%
% Inputs:
% varargin - variable number of arguments (between 1 and 4)
% --- varargin has the form: {object_type, opts, uinput, n_pde}
% Required:
% 1) varargin(1): data structure of the problem: PDE, DDE or PIE    
% --- object_type: mandatory input (PDE, DDE, or PIE structure)
% ----- PIE structure of the problem specifies PI operators, 
%       {T,Tu,Tw, A, Bi, Ci, Dij} as fields
% ----- if varargin(1) is PDE or DDE, then opts and uinput are optional 
%       (n_pde not used)
% ----- if varargin(1) is PIE, then opts/uinput/n_pde are mandatory

% --- varargin(2): opts - options for simulation parameters. 
%     If not provided, will be set to default values
% ------ opts is a structure with fields: 
%        {'N','tf','intScheme','Norder','dt','plot','type'}
% ------ N: order of polynomials used in spectral decomposition 
% ------ tf: final time of simulation
% ------ intScheme: time-integration scheme (backward difference scheme if 
%        1, symbolic integration of analytical solution if 2)
% ------ Norder: order of time-integration scheme (used only when 
%        intScheme=1, discretization order of backward difference scheme)
% ------ dt: time steps used in numerical time-integration
% ------ plot: returns plot if 1, no plot if 0
% ------ type: PDE, DDE, or PIE

% --- varargin(3): uinput - user-defined boundary inputs, forcing and 
%     initial conditions. If not provided, will be set to default values
% ------ uinput is a structure with fields: {'ic','w','ifexact','exact'}
% ------ ic: initial condition defined as a matlab symbolic object in 'sx' 
%        for each state in object_type separated as ic.ode and ic.pde
% ------ w: external disturbance defined as matlab symbolic object in 'st' 
%        and 'sx'
% ------ ifexact: if exact solution is known set to 1, else 0 (optional)
% ------ exact: only if uinput.ifexact=1, then define exact solution as a 
%        matlab symbolic object in 'sx' and 'st'

% --- varargin(4): n_pde (only used when varargin(1) is PIE), is a vector 
%     of integers specifying the number of differentiable states based on
%     index location, n_pde(i) is number of states that are (i-1) times
%     differentiable in space
%     For example, [1,2,3] stands for (1) continuous state, (2) continuously 
%     differentiable, and (3) twice continuously differentiable states 

% Outputs: 
% 1) solution 
% solution is a strucutre with the following fields
% --- solution.tf - actual final time of the solution
% --- solution.final.pde - pde (distributed state) solution at a final time : matrix of
%     dimension (N+1) x ns, ns=n0+n1+n2
% --- solution.final.ode - ode solution at a final time : array of dimensions
%     nx x 1

% IF OPTS.INTSCHEME=1 (BDF) OPTION (there are additional outputs)
% --- solution.timedep.ode - array of size nx x Nsteps - time-dependent solution of nx ODE states
% --- solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% --- solution of ns PDE (distributed) states of the primary PDE system
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5/26/2021
% Changing inputs and outputs to allow PIEs to be sent as inputs directly -
% SS 9/20/2021 - moved to PIESIM_options_check

function [solution, grid]=PIESIM(varargin)

% Check if options and uinput fields are defined
[opts, uinput]=PIESIM_options_check(varargin);

%--------------------------------------------------------------------
% Input check and conversion to PIE (if data structure is not PIE)
%--------------------------------------------------------------------

    if strcmp(opts.type,'PDE')
        PDE=varargin{1};
        % Check the inputs and define problem size
        [PDE, uinput, psize]=PIESIM_input_check(PDE, uinput, opts);
        
        % Convert PDE to PIE
        PIE=convert_PIETOOLS_PDE_batch(PDE);

    elseif strcmp(opts.type,'DDE')
        DDE=varargin{1};

        % Convert DDE to PIE
        DDE=initialize_PIETOOLS_DDE(DDE); % error checking and preprocessing
        DDF=minimize_PIETOOLS_DDE2DDF(DDE);
        PIE=convert_PIETOOLS_DDF2PIE(DDF);

        % Check the inputs and define problem size
        [PIE,uinput,psize]=PIESIM_input_check(PIE,uinput,opts);
   
    elseif strcmp(opts.type,'PIE')
        PIE=varargin{1};
        % Check the inputs and define problem size
        [PIE,uinput,psize]=PIESIM_input_check(PIE,uinput,opts);
    end  

    % Rescale all the PIE operators to be in domain [-1,1]
    PIE = rescalePIE(PIE,[-1,1]);
   
%-------------------------------------------------------------------------
% Set up initial conditions and boundary inputs for PIE from user-defined
% functions to the format required for PIE solution
%-------------------------------------------------------------------------

uinput=PIESIM_initial_setup(uinput,psize.n,opts.type);


% Setup a spatial discretization of the initial conditions and the PIE 
% operators with the Chebyshev Galerkin method

[Dop, coeff, grid]=PIESIM_discretize_all(PIE, uinput, psize);

% Solving in time for Chebyshev coefficients

disp('Setup completed: integrating in time');

solcoeff=PIESIM_time_integrate(psize, opts, uinput, coeff, Dop);

% Transform solution to physical space

solution=PIESIM_transform_to_solution(psize, PIE.Tu, PIE.Tw, Dop.Mcheb_nonsquare, uinput, grid, solcoeff, opts);


% Output and plot results

PIESIM_plot_solution(solution, psize, uinput,grid, opts)

end


