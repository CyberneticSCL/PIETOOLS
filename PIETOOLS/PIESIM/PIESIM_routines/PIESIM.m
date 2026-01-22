%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM.m     PIETOOLS 2025
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
% --- varargin has the form: {object_type, opts, uinput, ndiff}
% Required:
% 1) varargin(1): data structure of the problem: PDE, DDE or PIE    
% --- object_type: mandatory input (PDE, DDE, or PIE structure)
% ----- PIE structure of the problem specifies PI operators, 
%       {T,Tu,Tw, A, Bi, Ci, Dij} as fields
% ----- For any input structyre, opts, uinput and ndiff are
%       optional fields (defalut values will be used if not provided,
%       consult the manual)

% Optional:
% --- varargin(2): opts - options for simulation parameters. 
%     If not provided, will be set to default values
% ------ opts is a structure with fields: 
%        {'N','tf','intScheme','Norder','dt','plot','ploteig'}
% ------ N: order of polynomials used in spectral decomposition 
% ------ tf: final time of simulation
% ------ intScheme: time integration scheme (backward difference scheme if 
%        1, symbolic integration of analytical solution if 2)
% ------ Norder: order of time integration scheme (used only when 
%        intScheme=1, discretization order of backward difference scheme)
% ------ dt: time steps used in numerical time integration
% ------ plot: - option for plotting solution - returns solution plots if 'yes', no plots if 'no'
% ------ ploteig: option for plotting discrete eigenvalues of a discretized temporal propagator - returns eigenvalue plot if 'yes', no plot if 'no'
% ------ ifexact:  indicates whether a comparison of numerical solution with exact solution will be performed (if true) or not (if false). 
%        If uinput.ifexact=true and opts.plot=`yes', the exact solution will be plotted against numerical solution 
%------ dist: an optional flag that can take a value of `constant', `sin', or `sinc' - specifies a form of disturbance signal if disturbance is not defined by user. 
%------ control: an optional flag that can take a value of `constant', `sin', or `sinc' - specifies a form of control input signal if control input is not defined by user. 

% --- varargin(3): uinput - user-defined boundary inputs, forcing and 
%     initial conditions. If not provided, will be set to default values
% ------ uinput is a structure with fields: {'ic','w','ifexact','exact'}
% ------ ic: initial conditions defined as a MATLAB symbolic object. ic consists of the following sub-fields:
%        ic.pde - initial conditions for the PDE states defined as symbolic functions of 'sx' (1D), or 'sx', 'sy' (2D) 
%        ic.pie - initial conditions for the PIE states defined as symbolic functions of 'sx' (1D), or 'sx', 'sy' (2D) 
%        ic.ode - initial conditions for the ODE states defined as scalars
% ------ w: external disturbances defined as a MATLAB symbolic object in 'st','sx' (1D), or 'st','sx','sy' (2D)
% ------ exact: only if uinput.ifexact=true, then define exact solution as a 
%        MATLAB symbolic object in 'st','sx' (1D), or 'st','sx','sy' (2D)

% --- varargin(4): ndiff (only used when varargin(1) is PIE), is a vector 
%     of integers specifying the number of differentiable states based on
%     the index location, such as ndiff(i) is the number of states that are (i-1) times
%     differentiable in space.
%     For example, [0,2,1] stands for (0) continuous, (2) continuously 
%     differentiable, and (1) twice continuously differentiable states. If
%     not specified, the value of of 0 will be assumed for all
%     states 

% Outputs:
% 1) solution 
% solution is a structure with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.pde - array of size (N+1) x ns, ns is the total number of PDE states (including 2D and 1D states in 2D problems)- PDE (distributed state) solution at a final time
% --- solution.final.ode - array of size no - ode solution at a final time 
% --- solution.final.observed - array of size noo  - final value of observed outputs
% --- solution.final.regulated - array of size nro  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of discrete time values at which the time-dependent solution is computed  
% --- solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% --- solution.timedep.ode - array of size no x Nsteps - time-dependent solution of no ODE states
%     solution of ns PDE (distributed) states of the primary PDE system
% --- solution.timedep.observed - array of size noo x Nsteps -
%     time-dependent value of observed outputs
% --- solution.timedep.regulated - array of size nro x Nsteps -
%     time-dependent value of regulated outputs
%
%  2) grid - field containing two sub-fields
%  --- grid.phys 
%      a) In 1D - grid.phys is an array of size (N+1) x 1 containing spatial coordinates of the physical grid for
%  the primary solution
 %     b) In 2D - grid.phys is an array of size (N+1) x 2, with the first column containing spatial coordinates of the physical grid for
%  the primary solution along x direction, and the second column containing spatial coordinates of the physical grid for
%  the primary solution along y direction
%  --- grid.comp - in both 1D and 2D - array of size (N+1) x 1 containing spatial coordinates of the computaitonal grid along a single axis (mapped into [-1,1] domain) for
%  the primary solution
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
% YP 6/16/2022 - changed arguments in transform to solution to allow for different 
% LHS of PDE (T) and map (T0) operators (between fundamental and primary states) -
% needed, for example, for Orr-Sommerfeld equation
% DJ 10/10/2022 - Update for new terms format and command line parser.
% Batch format and old terms format can still be used, in which case
% opts.type is PDE_t or PDE_b.
% VJ 11/16/2025 - added support for cylindrical coordinate examples
% YP 1/22/2026 - passing u_tab, w_tab to initial setup

function [solution, grid]=PIESIM(varargin)

% Check if options and uinput fields are defined
[opts, uinput]=PIESIM_options_check(varargin);

%--------------------------------------------------------------------
% Input check and conversion to PIE (if data structure is not PIE)
%--------------------------------------------------------------------


    if strcmp(opts.type,'PDE') || strcmp(opts.type,'PDE_t') || strcmp(opts.type,'PDE_b')
        PDE = varargin{1};
        % Check the inputs and define problem size
        [PDE, uinput, opts, psize] = PIESIM_input_check(PDE, uinput, opts);
        
        % Convert PDE to PIE
        PIE = convert_PIETOOLS_PDE(PDE,'pie');

    elseif strcmp(opts.type,'DDE')
        DDE=varargin{1};

        % Convert DDE to PIE
        DDE = initialize_PIETOOLS_DDE(DDE); % error checking and preprocessing
        DDF = convert_PIETOOLS_DDE(DDE,'ddf');
        PIE = convert_PIETOOLS_DDF(DDF,'pie');

        % Check the inputs and define problem size
        [PIE,uinput,opts,psize]=PIESIM_input_check(PIE,uinput,opts);
   
    elseif strcmp(opts.type,'PIE')
        PIE=varargin{1};
        % Check the inputs and define problem size
        [PIE,uinput,opts,psize]=PIESIM_input_check(PIE,uinput,opts);
    end  
    
    % Rescale all the PIE operators to be in domain [-1,1]^d
    if (psize.dim==1)

    % Weighted Formulation for cylindrical PDE examples
    if isfield(uinput, 'weight') 
        PIE=weightPIE(PIE,uinput.weight);
    end
    
    PIE = rescalePIE(PIE,[-1,1]);
    else
          if isfield(uinput, 'weight') 
            disp('Warning: weighted formulation in 2D is not supported. Weight will be ignored.');
          end
    PIE = rescalePIE_2D(PIE,[-1,1;-1,1]);
    end
   
%-------------------------------------------------------------------------
% Set up initial conditions and boundary inputs for PIE from user-defined
% functions to the format required for PIE solution
%-------------------------------------------------------------------------

% Append u_tab and w_tab to psize for future use
psize.w_tab=PIE.w_tab;
psize.u_tab=PIE.u_tab;
uinput=PIESIM_initial_setup(uinput,psize,opts.type);


% Setup a spatial discretization of the initial conditions and the PIE 
% operators with the Chebyshev Galerkin method


if (psize.dim==1)
[Dop, coeff, grid]=PIESIM_discretize_all(PIE, uinput, psize);
else
[Dop, coeff, grid]=PIESIM_discretize_all_2D(PIE, uinput, psize);
end
% Solving in time for Chebyshev coefficients


disp('Setup completed: integrating in time');

solcoeff=PIESIM_time_integrate(psize, opts, uinput, coeff, Dop);



% Transform solution to physical space

if (psize.dim==1)
solution=PIESIM_transform_to_solution(psize, PIE, Dop, uinput, grid, solcoeff, opts);
else
solution=PIESIM_transform_to_solution_2D(psize, PIE, Dop, uinput, grid, solcoeff, opts);
end


% Output and plot results

if (psize.dim==1)
PIESIM_plot_solution(solution, psize, uinput,grid, opts);
else
PIESIM_plot_solution_2D(solution, psize, uinput,grid, opts);
end

end
