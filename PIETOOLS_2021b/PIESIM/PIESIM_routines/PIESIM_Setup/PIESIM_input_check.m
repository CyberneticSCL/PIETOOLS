%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_input_check.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [structure, uinput, psize]=PIESIM_input_check(varargin)
% Check if necessary inputs are defined to start the simulations
% NOTE: All other variables will be checked in PIETOOLS converter
% Inputs:
% varargin - variable number of arguments (between 1 and 4)
% Required:
% 1) varargin(1): data structure of the proglem: PDE, DDE or PIE
% PIE structure of the problem specifies PI operators, T,Tu,Tw, A, Bi, Ci, Dij as fields
% if varargin(1) is PDE or DDE, the rest of the inputs are optional
% if varargin(1) is PIE, the rest of the inputs are requires
% 2) varargin(2): opts - options for simulation parameters. If empty or incomplete, will be
% set to default values
% 3) varargin(3): uinput - user-defined boundary inputs, forcing and initial
% conditions. If empty or incomplete, will be set to default values
% Not used for PDE/DDE, required for PIE
% 4) varargin(4): n_pde - number of states with increasing differentiability, for example
% [1,2,3] stands for (1) continuous state, (2) continuously differentiable,
% and (3) twice continuously differentiable states - only used it data structure is PIE 

% Outputs: 
% 1) structure - PDE, DDE or PIE - updated data structure with variables properly defined (if
% previously undefined). 
% 2) uinput - updated user's input structure with variables properly defined (if
% previously undefined). 
% 3) psize - size of the problem. Includes nu, nw, nx, nf, N and n
% (corresponding to n_PDE)
% All properly defined variables are uchanged.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_1_2021

syms st sx;
structure=varargin{1};
uinput=varargin{2};
opts=varargin{3};

%-----------------------------------------------------------
% Input check for PDEs
%-----------------------------------------------------------

if (opts.type=='PDE')
    PDE=structure;

  % Checking of the PDE inputs begins


 if isfield(PDE,'n')

% Input is in terms format
 psize.n=PDE.n.n_pde;
 ns=length(PDE.n.n_pde)
 else

% Input is in batch format

 if ~isfield(PDE,'n0')
 disp('Warning: PDE.n0 is not defined. Defaulting to zero');
 PDE.n0=0;
 end
 if ~isfield(PDE,'n1')
 disp('Warning: PDE.n1 is not defined. Defaulting to zero');
 PDE.n1=0;
 end
 if ~isfield(PDE,'n2')
 disp('Warning: PDE.n2 is not defined. Defaulting to zero');
 PDE.n2=0;
 end

 psize.n=[PDE.n0 PDE.n1 PDE.n2];
 ns=PDE.n0+PDE.n1+PDE.n2;
 end

 if ~isfield(PDE,'nw')
 disp('Warning: PDE.nw is not defined. Defaulting to zero');
 PDE.nw=0;
 end
 if ~isfield(PDE,'nu')
 disp('Warning: PDE.nu is not defined. Defaulting to zero');
 PDE.nu=0;
 end
 if ~isfield(PDE,'nx')
 disp('Warning: PDE.nx is not defined. Defaulting to zero');
 PDE.nx=0;
 end
 
 if ~isfield(uinput,'ic')
 disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
 disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
 uinput.ic.PDE(1:ns)=sin(pi*sx);
 uinput.ic.ODE(1:PDE.nx)=0;
 end
 
 if ~isfield(uinput.ic,'PDE')
 disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
 uinput.ic.PDE(1:ns)=sin(pi*sx);
 elseif ~isempty(symvar(uinput.ic.PDE)) && any(~ismember(symvar(uinput.ic.PDE),{'sx'}))
        error('Initial conditions for PDE must be symbolic expressions in sx');
 end
 
 if (PDE.nx>0)
 if ~isfield(uinput.ic,'ODE')
 disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
 uinput.ic.ODE(1:PDE.nx)=0;
 elseif any(~isdouble(uinput.ic.ODE))
    error('Initial conditions for ODE must be scalar constants');
 end
 end
 
 if(size(uinput.ic.PDE,2)<ns)
 disp('Warning: Number of initial conditions on PDE states is less than the number of PDE states');
 disp('Defalting the rest to zero');
 uinput.ic.PDE(1,size(uinput.ic.PDE,2)+1:ns)=0;
 end
 
 if(size(uinput.ic.PDE,2)>ns)
 disp('Warning: Number of initial conditions on PDE states is greater than the number of PDE states');
 disp('Defaulting all initial conditions to zero');
 uinput.ic=rmfield(uinput.ic,'PDE');
 uinput.ic.PDE(1:ns)=0.;
 end
 
 if (PDE.nx>0) 
 if(size(uinput.ic.ODE,2)<PDE.nx)
 disp('Warning: Number of initial conditions on ODE states is less than the number of ODE states');
 disp('Defalting the rest to zero');
 uinput.ic.ODE(size(uinput.ic.PDE,2)+1:PDE.nx)=0;
 end
 
 if(size(uinput.ic.ODE,2)>PDE.nx)
 disp('Warning: Number of initial conditions on ODE states is greater than the number of ODE states');
 disp('Defaulting all ODE initial conditions to zero');
 clear uinput.ic.ODE;
 uinput.ic=rmfield(uinput.ic,'ODE');
 uinput.ic.ODE(1:PDE.nx)=0.;
 end
 end

 
 if (PDE.nu>0)
 if ~isfield(uinput,'u')
 disp('Warning: nu is greater than zero, but user-defiened u inputs are not provided. Defaulting PDE.nu to zero.');
 PDE.nu=0;
 end
 end
 
 if (PDE.nu>0)
 if (size(uinput.u,2)<PDE.nu)
 disp('Warning: Number of provided u inputs is less than nu.');    
 disp('Defalting the rest of u inputs and their time derivatives to zero');
 uinput.u(size(uinput.u,2)+1:PDE.nu)=0;
 uinput.udot(size(uinput.u,2)+1:PDE.nu)=0;
 end
 
 if (size(uinput.u,2)>PDE.nu)
 disp('Warning: Number of provided u inputs is  greater than nu.');    
 disp('Defalting PDE.nu to zero');
 PDE.nu=0;
 end
 end
 
 % Check disturbance setup
 if (PDE.nw>0)
 if ~isfield(uinput,'w')
 disp('Warning: nw is greater than zero, but user-defiened w inputs are not provided. Defaulting PDE.nw to zero.');
 PDE.nw=0;
 elseif ~isempty(symvar(uinput.w)) && any(~ismember(symvar(uinput.w),{'st'}))
        error('Disturbance inputs must be symbolic expressions in st');
 else
 if (size(uinput.w,2)<PDE.nw)
 disp('Warning: Number of provided w inputs is less than nw');    
 disp('Defalting the rest of w inputs and their time derivatives to zero');
 uinput.w(size(uinput.w,2)+1:PDE.nw)=0;
 uinput.wdot(size(uinput.w,2)+1:PDE.nw)=0;
 end
 if (size(uinput.w,2)>PDE.nw)
 disp('Warning: Number of provided w inputs is greater than nw');    
 disp('Defalting PDE.nw to zero');
 PDE.nw=0;
 end
 end

 if isfield(PDE,'Bw')
 nc_Bw=size(PDE.Bw,2);
 if (nc_Bw<PDE.nw)
 disp('Warning: Dimension of Bw matrix is incorrect. Initializing the missing columns to zero');  
 PDE.Bw(:,nc_Bw+1:PDE.nw)=0;
 end
 end


 if isfield(PDE,'B21_nonpol')   
 uinput.B21_nonpol=PDE.B21_nonpol;
 end


 end
 
 
 if ~isfield(uinput,'ifexact')   
 disp('Watning: uinput.ifexact is not specified. Defaulting  to false');
 uinput.ifexact=false;
 end
 
 if(uinput.ifexact)
 if ~isfield(uinput,'exact')
 disp('Warning: exact solution is not provided. Defaulting uinput.ifexact to false');
     uinput.ifexact=false;
 end 
 end
 
 if(uinput.ifexact)
 if (size(uinput.exact,2)~=ns)
 disp('Warning: number of exact solutions provided does not match the number of PDE states');
 disp('Defaulting uinput.ifexact to false');
     uinput.ifexact=false;
 end 
 end
 
 if ~isfield(PDE,'dom')
 disp('PDE domain is not defined. Defaulting to [-1, 1]');
 PDE.dom=[-1 1];
 end
 
 if (PDE.dom(1)==PDE.dom(2))
 disp('Warning: left and right ends of the domain are the same. Defaulting domain to [-1, 1]');
 PDE.dom=[-1 1];
 end

 uinput.a=PDE.dom(1);
 uinput.b=PDE.dom(2);

  % Checking of the PDE inputs ends
 structure=PDE;

 % Define problem size for discretization
        
  psize.nu=PDE.nu;
  psize.nw=PDE.nw;
  psize.nx=PDE.nx;
  psize.N=opts.N;
           
%-----------------------------------------------------------
% Input check for DDEs
%-----------------------------------------------------------
elseif (opts.type=='DDE')
        PIE=structure;

% Define problem size for discretization
        psize.nu=PIE.Tu.dim(1,2);
        psize.nw=PIE.Tw.dim(1,2);
        psize.nx=PIE.T.dim(1,1);
        ns=PIE.T.dim(2,1);
        psize.N=opts.N;
        psize.n=[0 ns];

  % Checking of the DDE inputs begins

        uinput.a=-1;
        uinput.b=0;
        uinput.ifexact=false;
        
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
            disp('Initial conditions on history are not defined. Defaulting to 1');
            uinput.ic.PDE(1:ns)=1;
        else
            if ~isfield(uinput.ic,'ODE')
                disp('Initial conditions on ODE states are not defined. Defaulting to 1');
                uinput.ic.ODE(1:psize.nx)=1;
            end
            if ~isfield(uinput.ic,'PDE')
                disp('Initial conditions on history are not defined. Defaulting to 1');
                uinput.ic.PDE(1:ns)=1;
            end
            
            if size(uinput.ic.ODE,1)>1
                uinput.ic.ODE=uinput.ic.ODE';
            end
            if max(size(uinput.ic.ODE))<psize.nx
                disp('Not enought number of initial conditions is defined. Defaulting the rest to 1');
                uinput.ic.ODE(max(size(uinput.ic.ODE)):psize.nx)=1;
            elseif max(size(uinput.ic.ODE))>psize.nx
                disp('Too many initial conditions are defined. Ignoring extra conditions');
                new_ic(1:psize.nx)=uinput.ic.ODE(1:psize.nx);
                uinput.ic.ODE=new_ic;
            end
            
            
        end

  % Checking of the DDE inputs ends

%-----------------------------------------------------------
% Input check for PIEs
%-----------------------------------------------------------

elseif (opts.type=='PIE')
    PIE=structure;

    % Define problem size for discretization    
    psize.nu=PIE.Tu.dim(1,2);
    psize.nw=PIE.Tw.dim(1,2);
    psize.nx=PIE.T.dim(1,1);
    psize.N=opts.N;
    psize.n=opts.piesize;
    
    if isfield(PIE,'dom')
        uinput.a=PIE.dom(1);
        uinput.b=PIE.dom(2);
    end
    
    
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
    
    if ~isfield(uinput,'w')
      uinput.w(1:psize.nw) = 0;
      disp('Warning: disturbances are not defined. Defaulting to zero');
    end

    if ~isfield(uinput,'u')
      uinput.u(1:psize.nu) = 0;
      disp('Warning: control inputs are not defined. Defaulting to zero');
    end
 
    if ~isfield(uinput,'ic')
        uinput.ic.ODE(1:PIE.T.dim(1,1))=1;
        uinput.ic.PDE(1:PIE.T.dim(2,1)) = sin(pi*sx);
        disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
        disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
    end
   

end
 
 
 

 
 
 
  



