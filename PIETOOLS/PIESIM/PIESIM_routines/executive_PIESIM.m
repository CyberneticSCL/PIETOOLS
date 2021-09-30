%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% executive_PIESIM.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please, contact Y. Peet at ypeet@asu.edu for support

% This rotuine is the executive routine for PIESIM. It performs the following functions:
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

% Inputs:
% Inputs can be specified in the following two ways:
% Input set one -
% structure - data structure of the problem (PDE or DDE)
% opts - options for simulation parameters
% uinput - user-defined boundary inputs
%
% Input set two -
% PIE - PIE structure of the problem, with PI operators, T,Tu,Tw, A, Bi, Ci, Dij as fields
% opts - options for temporal scheme parameters
% uinput - user-defined boundary inputs
% n_pde - number of states with increasing differentiability, for example
% [1,2,3] stands for (1) continuous state, (2) continuously differentiable,
% and (3) twice continuously differentiable states
%
% Outputs: solution. This field contains
% AVAILABLE FOR ALL OPTS.INTSCHEME OPTIONS
% solution.tf - actual final time of the solution
% solution.final.pde - pde (distributed state) solution at a final time : matrix of
% dimension (N+1) x ns, ns=n0+n1+n2
% solution.final.ode - ode solution at a final time : array of dimensions
% nx x 1

% AVAILABLE ONLY FOR OPTS.INTSCHEME=1 (BDF) OPTION
% solution.timedep.ode - array of size nx x Nsteps - time-dependent solution of nx ODE states
% solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% solution of ns PDE (distributed) states of the primary PDE system
% solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_26_2021
% Changing inputs and outputs to allow PIEs to be sent as inputs directyl -
% SS 9/20/2021

function solution=executive_PIESIM(varargin)

% required fields for options and uinputs
fields_opts = {'N','tf','intScheme','Norder','dt','Nsteps','plot'};
default_opts = {8, 10, 1, 2, 0.01, 1000, 'no'};

if nargin==1
    opts = struct();
    uinput = struct();
elseif nargin==2
    opts = varargin{2};
    uinput = struct();
else
    opts = varargin{2};
    uinput = varargin{3};
end

% check if all options are defined, if not define them
for i=1: length(fields_opts)
    if ~isfield(opts,fields_opts{i})
        opts.(fields_opts{i}) = default_opts{i};
    end
end

if isNOTPIE(varargin{1}) % if input type is PDE or DDE then convert to PIE first
    syms st;
    
    structure = varargin{1};
    
    if isfield(structure,'tau')
        disp('Solving DDE problem');
        DDE=structure;
        structure.type='DDE';
        opts.type='DDE';
    else
        disp('Solving PDE problem');
        PDE=structure;
        structure.type='PDE';
        opts.type='PDE';
    end
    
    
    % Check data structure for completeness
    if exist('PDE','var')
        [PDE, uinput]=PIESIM_input_check(PDE, uinput);
        
        % Transform PDE from the user-defined PDE domain x\in [a,b] to the computational domain xi\in[-1,1]
        
        % Place boundaries of the PDE domain into the variables "uinput.a" and "uinput.b"
        uinput.a=PDE.dom(1);
        uinput.b=PDE.dom(2);
        PDE=PIESIM_domain_transform(PDE);
        
        % Convert PDE to PIE
        
        PIE=convert_PIETOOLS_PDE_batch(PDE);
        % Define problem size for discretization
        
        psize.nu=PDE.nu;
        psize.nw=PDE.nw;
        psize.nx=PDE.nx;
        psize.nf=PDE.nf;
        psize.N=opts.N;
        psize.n=[PDE.n0 PDE.n1 PDE.n2];
    elseif exist('DDE','var')
        DDE=structure;
        DDE=initialize_PIETOOLS_DDE(DDE); % error checking and preprocessing
        DDF=minimize_PIETOOLS_DDE2DDF(DDE);
        PIE=convert_PIETOOLS_DDF2PIE(DDF);
        psize.nu=PIE.Tu.dim(1,2);
        psize.nw=PIE.Tw.dim(1,2);
        psize.nx=PIE.T.dim(1,1);
        psize.nf=0;
        ns=PIE.T.dim(2,1);
        psize.N=opts.N;
        psize.n=ns;
        uinput.a=-1;
        uinput.b=0;
        uinput.ifexact=false;
        % All this should go into subroutine
        
        if (psize.nw>0)
            if ~isfield(uinput,'w')
                disp('Disturbances are not defined. Defaulting to a user-specified signal type');
                if ~isfield(opts,'dist')
                    disp('Signal type for disturbances is not defined. Defaulting to a sinusoidal type');
                    uinput.w(1:psize.nw)=sin(st);
                else
                    switch opts.dist
                        case 'constant'
                            uinput.w(1:psize.nw)=1+0*st;
                        case 'sin'
                            uinput.w(1:psize.nw)=sin(st);
                        case 'sinc'
                            uinput.w(1:psize.nw)=sinc(st);
                        otherwise
                            disp('Signal type for disturbances is not defined. Defaulting to a sinusoidal type');
                            uinput.w(1:psize.nw)=sin(st);
                    end
                end
            end
        end
        if (psize.nu>0)
            if ~isfield(uinput,'u')
                disp('Control inputs are not defined. Defaulting to zero');
                uinput.u(1:psize.nu)=0;
            end
        end
        if ~isfield(uinput,'ic')
            disp('Initial conditions on ODE states are not defined. Defaulting to 1');
            uinput.ic.ODE(1:psize.nx)=1;
            disp('Initial conditions on history are not defined. Defaulting to linear');
            uinput.ic.PDE(1:ns)=1;
        else
            if ~isfield(uinput.ic,'ODE')
                disp('Initial conditions on ODE states are not defined. Defaulting to 1');
                uinput.ic.ODE(1:psize.nx)=1;
            end
            if ~isfield(uinput.ic,'PDE')
                disp('Initial conditions on history are not defined. Defaulting to linear');
                uinput.ic.PDE(1:ns)=1;
            end
            
            if size(uinput.ic.ODE,1)>1
                uinput.ic.ODE=uinput.ic.ODE';
            end
            if max(size(uinput.ic.ODE))<psize.nx
                disp('Not enought number of initial conditions is defined. Defaulting the rest to 1');
                uinput.ic.ODE(max(size(uinput.ic.ODE)):psize.nx)=1;
            elseif max(size(uinput.ic.ODE))>psize.nx
                disp('Too many initial conditions is defined. Ignoring extra conditions');
                new_ic(1:psize.nx)=uinput.ic.ODE(1:psize.nx);
                uinput.ic.ODE=new_ic;
            end
            
            
        end
        % Rescale all the PIE operators to be in domain [-1,1]
        PIE = rescalePIE(PIE,[-1,1]);
    end
    
    
elseif nargin==4
    PIE = varargin{1};
    
    % Define problem size for discretization
    
    psize.nu=PIE.Tu.dim(1,2);
    psize.nw=PIE.Tw.dim(1,2);
    psize.nx=PIE.T.dim(1,1);
    psize.nf= 0; % assuming no forcing, NEED to fix later
    psize.N=opts.N;
    psize.n=varargin{4};
    opts.type='PIE';
    opts.plot='no';
    
    % check if all uinputs are defined, if not define them
    if ~isfield(uinput,'a')
        uinput.a = 0;
    end
    if ~isfield(uinput,'b')
        uinput.b = 1;
    end
    if ~isfield(uinput,'ifexact')
        uinput.ifexact = false;
    end
    for i ={'w','u'}
        if ~isfield(uinput,i)
            uinput.(i{:}) = 0;
        end
    end
    if ~isfield(uinput,'ic')
        uinput.ic.ODE(1:PIE.T.dim(1,1))=1;
        uinput.ic.PDE(1:PIE.T.dim(2,1)) = 1;
    end
    
else
    error("If input object type is PIE, then number of differentiable states should be specified, for example 'executive_PIESIM(PIE,opts,uinput,n_pde)'");
end % end ifPDE




% Set up initial conditions and boundary inputs for PIE from user-defined
% functions to the format required for PIE solution

uinput=PIESIM_initial_setup(uinput);



% Setup a spatial discretization of the initial conditions and the PIE operators with the Chebyshev Galerkin method

[Dop, coeff, grid]=PIESIM_discretize_all(PIE, uinput, psize);

% Solving in time for Chebyshev coefficients

disp('Setup completed: integrating in time');
solcoeff=PIESIM_time_integrate(psize, opts, uinput, coeff, Dop);

% Transform solution to physical space

solution=PIESIM_transform_to_solution(psize, PIE.Tu, PIE.Tw, Dop.Mcheb_nonsquare, uinput, grid, solcoeff, opts);


% Output and plot results

PIESIM_plot_solution(solution, psize, uinput,grid, opts)

end

function logval = isNOTPIE(obj)
logval = ~isfield(obj,'T');
end

