%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_input_check.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [structure, uinput, opts, psize]=PIESIM_input_check(varargin)
% Check if necessary inputs are defined to start the simulations
% NOTE: All other variables will be checked in PIETOOLS converter
% Inputs:
% varargin - variable number of arguments (between 1 and 4)
% Required:
% 1) varargin{1}: data structure of the proglem: PDE, DDE or PIE
% PIE structure of the problem specifies PI operators, T,Tu,Tw, A, Bi, Ci, Dij as fields
% if varargin{1} is PDE or DDE, the rest of the inputs are optional
% if varargin{1} is PIE, ndiff input is required, while other inputs are
% optional
% 2) varargin{2}: opts - options for simulation parameters. If empty or incomplete, will be
% set to default values
% 3) varargin{3}: uinput - user-defined boundary inputs, forcing and initial
% conditions. If empty or incomplete, will be set to default values
% Not used for PDE/DDE, required for PIE
% 4) varargin{4}: ndiff - number of states with increasing differentiability, for example
% [1,2,3] stands for (1) continuous state, (2) continuously differentiable,
% and (3) twice continuously differentiable states - only used it data structure is PIE

% Outputs:
% 1) structure - PDE, DDE or PIE - updated data structure with variables properly defined (if
% previously undefined).
% 2) uinput - updated user's input structure with variables properly defined (if
% previously undefined).
% 3) opts: opts.ifexact will be updated to 'false' if exact solution is not
% provided
% 3) psize - size of the problem. 
% All properly defined variables are uchanged.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_1_2021
% YP 6/16/2022 - Separated inputs check for batch and terms format
% of PDE. Enhanced the input check for terms format. Terms format has now
% its own structure ('PDT'), while batch format has 'PDE' structure.
% Renamed uinput.B21_nonpol to uinput.Bpw_nonpol
% DJ, 10/10/2022 - Update to new terms format. The function can still
% process batch format and old terms format, but these functions have been
% added to separate subroutines.
%   Also set values for noo, nro in case of type=='DDE'
% DJ, 12/16/2024: Change default variables to (s1,s1_dum) in ODE case;
% DJ, 12/28/2024: Add spatial domain to "uinput" in 2D case;
% YP, 06/03/2025: Fixed sensing of observed and regulated outputs during coupling with PIE
% YP, 12/29/2025: Updated disturbance and control input check
% YP, 1/27/2026: Added support for infinite-dimensional outputs in 2D

syms st sx sy;
structure=varargin{1};
uinput=varargin{2};
opts=varargin{3};


if strcmp(opts.type,'PDE_b')
    %-----------------------------------------------------------
    % Input check for PDEs in batch format
    %-----------------------------------------------------------
    [structure, uinput, opts, psize] = PIESIM_input_check_PDE_batch(varargin{:});
    psize.dim = 1;

elseif strcmp(opts.type,'PDE_t')
    %-----------------------------------------------------------
    % Input check for PDEs in legacy terms format
    %-----------------------------------------------------------
    [structure, uinput, opts, psize] = PIESIM_input_check_PDE_terms_legacy(varargin{:});
    psize.dim = 1

elseif strcmp(opts.type,'PDE')
    %-----------------------------------------------------------
    % Input check for PDEs in terms format
    %-----------------------------------------------------------
    if isa(structure,'struct')
        structure = pde_struct(structure);  % Convert struct to pde_struct
    elseif isa(structure,'sys')
        structure = structure.params;   % Extract terms from sys object
    end

    PDE = initialize(structure);
    [PDE,x_order] = reorder_comps(PDE,'x'); % Reorder components in increasing order of differentiability
    x_order = vec_order2elem_order(x_order,PDE.x_tab(:,2)); % Adjust order to account for vector-valued state components;

    psize = struct();
    psize.dim = PDE.dim;
    % Checking of the PDE inputs begins


    if PDE.dim==0
        % If the system is finite-dimensional, augment to 1D
        pvar s1 s1_dum                                                      % DJ, 12/16/2024
        vars = [s1,s1_dum];
        dom = [-1,1];

        % Add a temporary state variable that exists on the 1D domain.
        ncomps = numel(PDE.x);
        PDE.x{ncomps+1}.dom = dom;
        PDE.x{ncomps+1}.vars = vars;
        PDE.x{ncomps+1}.term{1}.x = ncomps+1;
        % Initialize the augmented system, and get rid of the temporary state.
        PDE = initialize_PIETOOLS_PDE(psize,true);
        PDE.x = PDE.x((1:ncomps)');
        PDE.x_tab = PDE.x_tab((1:ncomps)',:);
    end

% % Establish size of ODE and PDE state components
x_tab = PDE.x_tab;
% Determine the number of ODE state variables
psize.no = sum(x_tab(~any(x_tab(:,3:2+PDE.dim),2),2));
if PDE.dim==1
    % Determine the number of PDE state variables differentiable up to each
    % order
    nmax = max(x_tab(:,4));
    psize.n = zeros(1,[nmax+1]);
    for k=0:nmax
        psize.n(k+1) = sum(x_tab(x_tab(:,3) & x_tab(:,4)==k,2));
    end
elseif PDE.dim==2
    % Determine which state variables vary in either x, y, or both
    is1D_x = all(x_tab(:,3:2+PDE.dim)==[1,0],2);
    is1D_y = all(x_tab(:,3:2+PDE.dim)==[0,1],2);
    is2D = all(x_tab(:,3:2+PDE.dim),2);
    if ~any(is1D_x)
        psize.nx = zeros(1,0);
    else
        % Determine the maximal order of differentiability along x-direction
        nmax_x = max(x_tab(is1D_x,5));
        
        % Determine the number of state variables differentiable up to each
        % order
        psize.nx = zeros(1,nmax_x+1);
        for k=0:nmax_x
            psize.nx(k+1) = sum(x_tab(is1D_x & x_tab(:,5)==k,2));
        end
    end
    if ~any(is1D_y)
        psize.ny = zeros(1,0);
    else
        % Determine the maximal order of differentiability along y-direction
        nmax_y = max(x_tab(is1D_y,6));
        
        % Determine the number of state variables differentiable up to each
        % order
        psize.ny = zeros(1,nmax_y+1);
        for k=0:nmax_y
            psize.ny(k+1) = sum(x_tab(is1D_y & x_tab(:,6)==k,2));
        end
    end
    
    if ~any(is2D)
        psize.n = zeros([0,0]);        
    else
        % Determine the maximal order of differentiability along each direction
        nmax_2D = max(x_tab(is2D,[5,6]),[],1);
        
        % Determine the number of state variables differentiable up to each
        % order
        psize.n = zeros(nmax_2D+1);
        for k=1:numel(psize.n)
            [i,j] = ind2sub(nmax_2D+1,k);
            psize.n(k) = sum(x_tab(is2D & all(x_tab(:,[5,6])==[i-1,j-1],2),2));
        end
    end
end


% Convert psize.n back to array for now

if PDE.dim==2
    nmax=max(size(psize.n));
% 
    array=zeros(1,nmax);  
% 
    for i=1:size(psize.n,1)
        for j=1:size(psize.n,2)
            ind=max(i,j);
            array(ind)=array(ind)+psize.n(i,j);
        end
    end
% 
    psize.n=array;
     end % if  PDE.dim==2

    % Split vector-valued state components into separate scalar-valued ones
    psize.x_tab = zeros(sum(x_tab(:,2)),size(x_tab,2));
    r_strt = 0;
    for i=1:size(x_tab,1)
        nr = x_tab(i,2);    % Check the size of state component i
        psize.x_tab(r_strt+(1:nr),:) = repmat(x_tab(i,:),nr,1);
        r_strt = r_strt+nr;
    end
    psize.dim=PDE.dim;

%    x_order_PDE = x_order(psize.no+1:end);
    if PDE.dim==1   
    ns = sum(psize.n,'all');
    else
        nsx=sum(psize.nx,'all');
        nsy=sum(psize.ny,'all');
        ns2d=sum(psize.n,'all');
    ns = nsx+nsy+ns2d;
    end

    % Establish sizes of inputs and outputs
    psize.nw = sum(PDE.w_tab(:,2)); % number of disturbances
    psize.nu = sum(PDE.u_tab(:,2)); % number of control inputs

   % Define sizes for inputs and outputs
    % Determine which inpus and outputs vary in either x, y, or both

     if (PDE.dim==1)

    psize.nwx = sum(PDE.w_tab(:,3) == 1);
    psize.nwy = 0;
    psize.nw2 = 0;
    psize.nw0 = psize.nw-(psize.nwx+psize.nwy+psize.nw2);

    psize.nux = sum(PDE.u_tab(:,3) == 1);
    psize.nuy = 0;
    psize.nu2 = 0;
    psize.nu0 = psize.nu-(psize.nux+psize.nuy+psize.nu2);

    psize.nrox=sum(PDE.z_tab(:,3)); % number of infinite-dimensional regulated outputs
    psize.nro = sum(PDE.z_tab(:,2))-psize.nrox; % number of finite-dimensional regulated outputs
    psize.noox = sum(PDE.y_tab(:,3)); % number of infinite-dimensional observed outputs
    psize.noo = sum(PDE.y_tab(:,2))-psize.noox; % number of finite-dimensional observed outputs
    % % 
     else

    psize.nwx = sum(PDE.w_tab(:,3) == 1 & PDE.w_tab(:,2+PDE.dim) == 0);
    psize.nwy = sum(PDE.w_tab(:,3) == 0 & PDE.w_tab(:,2+PDE.dim) == 1);
    psize.nw2 = sum(PDE.w_tab(:,3) == 1 & PDE.w_tab(:,2+PDE.dim) == 1);
    psize.nw0 = psize.nw-(psize.nwx+psize.nwy+psize.nw2);

    psize.nux = sum(PDE.u_tab(:,3) == 1 & PDE.u_tab(:,2+PDE.dim) == 0);
    psize.nuy = sum(PDE.u_tab(:,3) == 0 & PDE.u_tab(:,2+PDE.dim) == 1);
    psize.nu2 = sum(PDE.u_tab(:,3) == 1 & PDE.u_tab(:,2+PDE.dim) == 1);
    psize.nu0 = psize.nu-(psize.nux+psize.nuy+psize.nu2);

    psize.nrox = sum(PDE.z_tab(:,3) == 1 & PDE.z_tab(:,2+PDE.dim) == 0);
    psize.nroy = sum(PDE.z_tab(:,3) == 0 & PDE.z_tab(:,2+PDE.dim) == 1);
    psize.nro2 = sum(PDE.z_tab(:,3) == 1 & PDE.z_tab(:,2+PDE.dim) == 1);
    psize.nro = sum(PDE.z_tab(:,2))-(psize.nrox+psize.nroy+psize.nro2);

    psize.noox = sum(PDE.y_tab(:,3) == 1 & PDE.y_tab(:,2+PDE.dim) == 0);
    psize.nooy = sum(PDE.y_tab(:,3) == 0 & PDE.y_tab(:,2+PDE.dim) == 1);
    psize.noo2 = sum(PDE.y_tab(:,3) == 1 & PDE.y_tab(:,2+PDE.dim) == 1);
    psize.noo = sum(PDE.y_tab(:,2))-(psize.noox+psize.nooy+psize.noo2);

  

     end

    % Define problem size for discretization
    psize.N = opts.N;

    if ~isfield(uinput,'ic')
        if psize.no>0
            disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
        end
        if ns>0
            disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
        end
        uinput.ic.ODE(1:psize.no)=0;
        if (PDE.dim==1)
        uinput.ic.PDE(1:ns)=sin(pi*sx);
        elseif (PDE.dim==2)
            uinput.ic.PDE(1:nsx)=sin(pi*sx);
            uinput.ic.PDE(nsx+1:nsx+nsy)=sin(pi*sy);
            uinput.ic.PDE(nsx+nsy+1:ns)=sin(pi*sx)*sin(pi*sy);
        end
    end

    ireorder=0;

    if isfield(uinput.ic,'x')
        if (isfield(uinput.ic,'ODE') || isfield(uinput.ic,'PDE'))
            disp('Initial conditions have been specified using both "ic.x" and "ic.ODE", "ic.PDE". Continuing with the initial conditions "ic.x".')
        end
        if ~isempty(symvar(sym(uinput.ic.x))) && any(~ismember(symvar(uinput.ic.x),{'sx';'sy'}))
            error('Initial conditions for PDE must be symbolic expressions in sx (and sy for 2D PDEs)');
        elseif size(uinput.ic.x,2)~=(psize.no + ns)
            error('Number of initial conditions should match the number of state variables')
        end
        if ~issorted(x_order)
        uinput.ic.x = uinput.ic.x(x_order); % Reorder initial conditions to match new ordering of state components.
        ireorder=1;
        disp('Initial conditions have been reordered to correspond to new ordering');
        end
        if hasSymType(uinput.ic.x(1:psize.no),'variable')
            error('Initial conditions for ODE must be scalar constants');
        else
            uinput.ic.ODE = double(uinput.ic.x(1:psize.no));
        end
        uinput.ic.PDE = uinput.ic.x(psize.no+1:end);
    end

       % Check initial conditions for ODE state variables
    if (psize.no>0)
        if ~isfield(uinput.ic,'ODE') || isempty(uinput.ic.ODE)
            disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
            uinput.ic.ODE(1:psize.no)=0;
        elseif any(~isdouble(uinput.ic.ODE))
            error('Initial conditions for ODE must be scalar constants');
        end
        if(size(uinput.ic.ODE,2)<psize.no)
            disp('Warning: Number of initial conditions on ODE states is less than the number of ODE states');
            disp('Defalting the rest to zero');
            uinput.ic.ODE(size(uinput.ic.ODE,2)+1:psize.no)=0;
        elseif(size(uinput.ic.ODE,2)>psize.no)
            disp('Warning: Number of initial conditions on ODE states is greater than the number of ODE states');
            disp('Defaulting all ODE initial conditions to zero');
            clear uinput.ic.ODE;
            uinput.ic=rmfield(uinput.ic,'ODE');
            uinput.ic.ODE(1:psize.no)=0;
        end
    end
    
    % Check initial conditions for PDE state variables

    if isfield(uinput.ic,'PDE')
    if(size(uinput.ic.PDE,2)<ns)
        disp('Warning: Number of initial conditions on PDE states is less than the number of PDE states');
        disp('Defaulting the rest to zero');
        uinput.ic.PDE(1,size(uinput.ic.PDE,2)+1:ns)=sym(0);
    elseif(size(uinput.ic.PDE,2)>ns)
        disp('Warning: Number of initial conditions on PDE states is greater than the number of PDE states');
        disp('Defaulting all initial conditions to zero');
        uinput.ic=rmfield(uinput.ic,'PDE');
        uinput.ic.PDE(1:ns)=sym(0);
    end
    % Reorder initial conditions only if it has bot been reordered already
        if ~issorted(x_order) & ~ireorder
            uinput.ic.x = sym(zeros(1,psize.no+ns));
            uinput.ic.x(1:psize.no)=sym(uinput.ic.ODE);
            uinput.ic.x(psize.no+1:psize.no+ns)=sym(uinput.ic.PDE); 
            uinput.ic.x = uinput.ic.x(x_order); % Reorder initial conditions to match new ordering of state components.
            uinput.ic.ODE = uinput.ic.x(1:psize.no);
            uinput.ic.PDE = uinput.ic.x(psize.no+1:end);
            disp('Initial conditions have been reordered to correspond to new ordering');
        end
    end

    if (PDE.dim==1)
    if ns>0
        if ~isfield(uinput.ic,'PDE')
            if isfield(uinput.ic,'PIE')
            disp('Initial conditions are provided as uinput.ic.PIE and not uinput.ic.PDE. They will be used for the PDE states as given.')
            uinput.ic.PDE=uinput.ic.PIE;
            else
            disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
            uinput.ic.PDE(1:ns)=sin(pi*sx);
            end
        elseif ~isempty(symvar(sym(uinput.ic.PDE))) && any(~ismember(symvar(uinput.ic.PDE),{'sx'}))
            error('Initial conditions for PDE must be symbolic expressions in sx');
        end
    end
    elseif(PDE.dim==2)
        if ns>0
        if ~isfield(uinput.ic,'PDE')
            disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
            uinput.ic.PDE(1:nsx)=sin(pi*sx);
            uinput.ic.PDE(nsx+1:nsx+nsy)=sin(pi*sy);
            uinput.ic.PDE(nsx+nsy+1:ns)=sin(pi*sx)*sin(pi*sy);
        end
        end
    end % PDE.dim

    if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.no+ns)
            disp('Warning: number of exact solutions provided does not match the number of PDE states');
            disp('Defaulting opts.ifexact to false');
            opts.ifexact=false;
        else
            if ~issorted(x_order)
        uinput.exact = uinput.exact(x_order); % Reorder exact solution to match new ordering of state components.
        disp('Exact solution has been reordered to correspond to new ordering');
            end
        end % ~isfield(uinput,'exact')
    end % opts.ifexact
   
    % Check disturbances and control inputs

    if isfield(opts,'dist')
        uinput.dist=opts.dist;
    end

     if isfield(opts,'control')
        uinput.control=opts.control;
    end

        uinput = PIESIM_disturbance_check(PDE,uinput,psize);
 
    if (PDE.dim<2)
    if (PDE.dom(1)==PDE.dom(2))
        disp('Warning: left and right ends of the domain are the same. Defaulting domain to [-1, 1]');
        PDE.dom=[-1 1];
    end
        uinput.a = PDE.dom(1);
        uinput.b = PDE.dom(2);
    else
        for k=1:2
        if (PDE.dom(k,1)==PDE.dom(k,2))
        disp('Warning: one of the domain segments is zero. Defaulting this segment to [-1, 1]');
        PDE.dom(k,:)=[-1 1];
        end
        end
        if ~isfield(uinput,'dom')                                           % DJ, 12/28/2024
            uinput.dom = PDE.dom;
        end
    end

    % Checking of the PDE inputs ends
    structure = PDE;


elseif (opts.type=='DDE')
    %-----------------------------------------------------------
    % Input check for DDEs
    %-----------------------------------------------------------
    PIE=structure;  % In DDE case, systems must be passed as PIE.

    % Define problem size for discretization
    psize.dim = 1;      % DDE is 1D
    psize.nu=size(PIE.Tu,2);
    psize.nw=size(PIE.Tw,2);
    psize.no=PIE.T.dim(1,1);
    ns = PIE.T.dim(2,1);
    psize.N=opts.N;
    psize.n=[0 ns];
    %------------------------
    % Size of finite-dimensional and infinite-dimensional outputs
    psize.nro = size(PIE.C1.Q1,1);
    psize.noo = size(PIE.C2.Q1,1); 
    psize.nrox = size(PIE.C1,1)-psize.nro;
    psize.noox = size(PIE.C2,1)-psize.noo;  

    psize.nu=size(PIE.Tu,2); % number of disturbances
    psize.nw=size(PIE.Tw,2); % number of control inputs

    % Determine which disturbances vary in either x, y, or both

    if ~isempty(PIE.w_tab)
    psize.nwx = sum(PIE.w_tab(:,3) == 1);
    else
    psize.nwx=0;
    end
    psize.nwy = 0;
    psize.nw2 = 0;
    psize.nw0 = psize.nw-(psize.nwx+psize.nwy+psize.nw2);

    if ~isempty(PIE.u_tab)
    psize.nux = sum(PIE.u_tab(:,3) == 1);
    else
    psize.nux=0;
    end
    psize.nuy = 0;
    psize.nu2 = 0;
    psize.nu0 = psize.nu-(psize.nux+psize.nuy+psize.nu2);
    % % 

    %------------------------

    % Checking of the DDE inputs begins

    uinput.a=-1;
    uinput.b=0;

    % Check disturbances and control inputs

    % Append opts.dist to unput.dist if exists

    if isfield(opts,'dist')
        uinput.dist=opts.dist;
    end

    if isfield(opts,'control')
        uinput.control=opts.control;
    end

    uinput = PIESIM_disturbance_check(PIE,uinput,psize);


    if ~isfield(uinput,'ic')
        if (psize.no>0)
        disp('Initial conditions on ODE states are not defined. Defaulting to 1');
        uinput.ic.ODE(1:psize.no)=1;
        end
        if ns>0
        disp('Initial conditions on history are not defined. Defaulting to 1');
        uinput.ic.DDE(1:ns)=1;
        end
    else
        if ~isfield(uinput.ic,'ODE')
            if psize.no>0
            disp('Initial conditions on ODE states are not defined. Defaulting to 1');
            uinput.ic.ODE(1:psize.no)=1;
            end
        end
        if ~isfield(uinput.ic,'DDE')
            if isfield(uinput.ic,'PDE') 
                disp('Initial conditions are provided as uinput.ic.PDE and not uinput.ic.DDE. They will be used for the DDE states as given.');
            uinput.ic.DDE=uinput.ic.PDE;
            else
            disp('Initial conditions on history are not defined. Defaulting to 1');
            uinput.ic.DDE(1:ns)=1;
            end

        end

        if(size(uinput.ic.DDE,2)<ns)
         disp('Warning: Number of initial conditions on DDE states is less than the number of DDE states');
         disp('Defalting the rest to zero');
         uinput.ic.DDE(1,size(uinput.ic.DDE,2)+1:ns)=0;
      end

if(size(uinput.ic.DDE,2)>ns)
    disp('Warning: Number of initial conditions on DDE states is greater than the number of DDE states');
    disp('Defaulting all initial conditions to zero');
    uinput.ic=rmfield(uinput.ic,'DDE');
    uinput.ic.DDE(1:ns)=0.;
end


        if size(uinput.ic.ODE,1)>1
            uinput.ic.ODE=uinput.ic.ODE';
        end
        if max(size(uinput.ic.ODE))<psize.no
            disp('Not enought number of ODE initial conditions is defined. Defaulting the rest to 1');
            uinput.ic.ODE(max(size(uinput.ic.ODE)):psize.no)=1;
        elseif max(size(uinput.ic.ODE))>psize.no
            disp('Too many ODE initial conditions are defined. Ignoring extra conditions');
            new_ic(1:psize.no)=uinput.ic.ODE(1:psize.no);
            uinput.ic.ODE=new_ic;
        end


    end

     if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.no+ns)
            disp('Warning: number of exact solutions provided does not match the number of PDE states');
            disp('Defaulting opts.ifexact to false');
            opts.ifexact=false;
        end % ~isfield(uinput,'exact')
    end % opts.ifexact

% Reassign uinput.ic.DDE to uinput.ic.PDE

uinput.ic.PDE=uinput.ic.DDE;

    % Checking of the DDE inputs ends

elseif (opts.type=='PIE')
    %-----------------------------------------------------------
    % Input check for PIEs
    %-----------------------------------------------------------
    PIE=structure;

    if isa(PIE,'pie_struct') || isfield(PIE,'dim')
        psize.dim = PIE.dim;
    elseif isa(PIE.T,'opvar2d')
        psize.dim = 2;
    else
        psize.dim = 1;
    end

    % Define problem size for discretization
    psize.no=PIE.T.dim(1,1);
    ns = PIE.T.dim(2,1);

% Find the number of finite-dimensional and infinite-dimensional outputs
% nro, noo counts finite-dimensional outputs
% nrox, noox counts infinite-dimensional outputs

    %------------------------
    psize.nro = size(PIE.C1.Q1,1);
    psize.noo = size(PIE.C2.Q1,1); 
    psize.nrox = size(PIE.C1,1)-psize.nro;
    psize.noox = size(PIE.C2,1)-psize.noo;  

    psize.nu=size(PIE.Tu,2); % number of disturbances
    psize.nw=size(PIE.Tw,2); % number of control inputs

    % Determine which disturbances vary in either x, y, or both

     if (psize.dim==1)
    
    if ~isempty(PIE.w_tab)     
    psize.nwx = sum(PIE.w_tab(:,3) == 1);
    else
    psize.nwx=0;
    end
    psize.nwy = 0;
    psize.nw2 = 0;
    psize.nw0 = psize.nw-(psize.nwx+psize.nwy+psize.nw2);

    if ~isempty(PIE.w_tab)   
    psize.nux = sum(PIE.u_tab(:,3) == 1);
    else
    psize.nux=0;
     end
    psize.nuy = 0;
    psize.nu2 = 0;
    psize.nu0 = psize.nu-(psize.nux+psize.nuy+psize.nu2);
    % % 
     else

    psize.nwx = sum(PIE.w_tab(:,3) == 1 & PIE.w_tab(:,2+psize.dim) == 0);
    psize.nwy = sum(PIE.w_tab(:,3) == 0 & PIE.w_tab(:,2+psize.dim) == 1);
    psize.nw2 = sum(PIE.w_tab(:,3) == 1 & PIE.w_tab(:,2+psize.dim) == 1);
    psize.nw0 = psize.nw-(psize.nwx+psize.nwy+psize.nw2);

    psize.nux = sum(PIE.u_tab(:,3) == 1 & PIE.u_tab(:,2+psize.dim) == 0);
    psize.nuy = sum(PIE.u_tab(:,3) == 0 & PIE.u_tab(:,2+psize.dim) == 1);
    psize.nu2 = sum(PIE.u_tab(:,3) == 1 & PIE.u_tab(:,2+psize.dim) == 1);
    psize.nu0 = psize.nu-(psize.nux+psize.nuy+psize.nu2);

     end

    %------------------------

    psize.N=opts.N;
    if ~isempty(opts.piesize)
    psize.n=opts.piesize;
    else
    psize.n=ns;
    end


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
  
     if isfield(opts,'dist')
        uinput.dist=opts.dist;
     end

     if isfield(opts,'control')
        uinput.control=opts.control;
     end  

    uinput = PIESIM_disturbance_check(PIE,uinput,psize);
    
    if ~isfield(uinput,'ic')
        if (psize.no>0)
        disp('Initial conditions on ODE states are not defined. Defaulting to 1');
        uinput.ic.ODE(1:psize.no)=1;
        end
        if (ns>0)
        disp('Initial conditions on PIE states  are not defined. Defaulting to 1');
        uinput.ic.PIE(1:ns)=1;
        end
    else
        if ~isfield(uinput.ic,'ODE')
            if (psize.no>0)
            disp('Initial conditions on ODE states are not defined. Defaulting to 1');
            uinput.ic.ODE(1:psize.no)=1;
            end
        end
        if ~isfield(uinput.ic,'PIE')
            if isfield(uinput.ic,'PDE')
            disp('Initial conditions are provided as uinput.ic.PDE and not uinput.ic.PIE. They will be used for the PIE states as given.');
            uinput.ic.PIE=uinput.ic.PDE;
        else
            disp('Initial conditions on PIE states are not defined. Defaulting to 1');
            uinput.ic.PIE(1:ns)=1;
        end
    end

        if isfield(uinput.ic,'ODE')
        if size(uinput.ic.ODE,1)>1
            uinput.ic.ODE=uinput.ic.ODE';
        end
        if max(size(uinput.ic.ODE))<psize.no
            disp('Not enought number of ODE initial conditions is defined. Defaulting the rest to 1');
            uinput.ic.ODE(max(size(uinput.ic.ODE)):psize.no)=1;
        elseif max(size(uinput.ic.ODE))>psize.no
            disp('Too many ODE initial conditions are defined. Ignoring extra conditions');
            new_ic(1:psize.no)=uinput.ic.ODE(1:psize.no);
            uinput.ic.ODE=new_ic;
        end
        end
    end

 if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.no+ns)
            disp('Warning: number of exact solutions provided does not match the number of PDE states');
            disp('Defaulting opts.ifexact to false');
            opts.ifexact=false;
        end % ~isfield(uinput,'exact')
    end % opts.ifexact

end

end


function [structure, uinput, opts, psize] = PIESIM_input_check_PDE_batch(varargin)
%-----------------------------------------------------------
%  Input check for PDEs in batch format - 1D cases only 
%  (infinite-dimensional inputs and outputs not supported)
% -----------------------------------------------------------

syms st sx;
structure=varargin{1};
uinput=varargin{2};
opts=varargin{3};

PDE=structure;
% Checking of the PDE inputs begins

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
ns = PDE.n0+PDE.n1+PDE.n2;

if ~isfield(PDE,'no')
    disp('Warning: PDE.no is not defined. Defaulting to zero');
    PDE.no=0;
end
if ~isfield(PDE,'nw')
    disp('Warning: PDE.nw is not defined. Defaulting to zero');
    PDE.nw=0;
end
if ~isfield(PDE,'nu')
    disp('Warning: PDE.nu is not defined. Defaulting to zero');
    PDE.nu=0;
end
if ~isfield(PDE,'nb')
    disp('Warning: PDE.nb is not defined. Defaulting to zero');
    PDE.nb=0;
end
if ~isfield(PDE,'nro')
    disp('Warning: PDE.nro is not defined. Defaulting to zero');
    PDE.nro=0;
end
if ~isfield(PDE,'noo')
    disp('Warning: PDE.noo is not defined. Defaulting to zero');
    PDE.noo=0;
end

if ~isfield(uinput,'ic')
    disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
    disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
    uinput.ic.PDE(1:ns)=sin(pi*sx);
    uinput.ic.ODE(1:PDE.no)=0;
end

if ~isfield(uinput.ic,'PDE')
    disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
    uinput.ic.PDE(1:ns)=sin(pi*sx);
elseif ~isempty(symvar(sym(uinput.ic.PDE))) && any(~ismember(symvar(uinput.ic.PDE),{'sx'}))
    error('Initial conditions for PDE must be symbolic expressions in sx');
end

if (PDE.no>0)
    if ~isfield(uinput.ic,'ODE')
        disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
        uinput.ic.ODE(1:PDE.no)=0;
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

if (PDE.no>0)
    if(size(uinput.ic.ODE,2)<PDE.no)
        disp('Warning: Number of initial conditions on ODE states is less than the number of ODE states');
        disp('Defalting the rest to zero');
        uinput.ic.ODE(size(uinput.ic.PDE,2)+1:PDE.no)=0;
    end

    if(size(uinput.ic.ODE,2)>PDE.no)
        disp('Warning: Number of initial conditions on ODE states is greater than the number of ODE states');
        disp('Defaulting all ODE initial conditions to zero');
        clear uinput.ic.ODE;
        uinput.ic=rmfield(uinput.ic,'ODE');
        uinput.ic.ODE(1:PDE.no)=0.;
    end
end

% Check control inputs setup
if (PDE.nu>0)
    if ~isfield(uinput,'u')
        disp('Warning: nu is greater than zero, but user-defiened u inputs are not provided. Defaulting PDE.nu to zero.');
        PDE.nu=0;
    else

    % Make sure uinput.u is a cell array
    if ~iscell(uinput.u)
        uinput.u=num2cell(uinput.u);
    end
        
        if (size(uinput.u,2)<PDE.nu)
            disp('Warning: Number of provided u inputs is less than nu');
            disp('Defalting the rest of u inputs to zero');
            uinput.u{size(uinput.u,2)+1:PDE.nu}=0;
        end
        if (size(uinput.u,2)>PDE.nu)
            disp('Warning: Number of provided u inputs is greater than nu');
            disp('Defalting PDE.nu to zero');
            PDE.nu=0;
        end
    end

end


% Check disturbance setup
if (PDE.nw>0)
    if ~isfield(uinput,'w')
        disp('Warning: nw is greater than zero, but user-defiened w inputs are not provided. Defaulting PDE.nw to zero.');
        PDE.nw=0;
    else

    % Make sure uinput.w is a cell array
    if ~iscell(uinput.w)
        uinput.w=num2cell(uinput.w);
    end
        
        if (size(uinput.w,2)<PDE.nw)
            disp('Warning: Number of provided w inputs is less than nw');
            disp('Defalting the rest of w inputs to zero');
            uinput.w(size(uinput.w,2)+1:PDE.nw)={0};
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
        uinput.Bpw_nonpol=PDE.B21_nonpol;
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


PDE.nx=PDE.no;
PDE.ny=PDE.noo;
PDE.nz=PDE.nro;

structure=PDE;

% Define problem size for discretization

psize.nu=PDE.nu; % number of control inputs
psize.nw=PDE.nw; % number of disturbances
psize.no=PDE.no; % number of ODE states
psize.nro=PDE.nro; % number of regulated outputs
psize.noo=PDE.noo; % number of observed outputs
psize.N=opts.N; % order of disceretization in space
psize.dim=1;


 if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.no+ns)
            disp('Warning: number of exact solutions provided does not match the number of PDE states');
            disp('Defaulting opts.ifexact to false');
            opts.ifexact=false;
        end % ~isfield(uinput,'exact')
    end % opts.ifexact

end


function [structure, uinput, opts, psize] = PIESIM_input_check_PDE_terms_legacy(varargin)
%-----------------------------------------------------------
%  Input check for PDEs in legacy terms format - only for 1D cases
%  (infinite-dimensional inputs and outputs not supported)
% -----------------------------------------------------------

syms st sx;
structure=varargin{1};
uinput=varargin{2};
opts=varargin{3};

psize=structure;
psize.dim=1;

% Define problem size for discretization

if ~isfield(psize.n,'nu')
    disp('Warning: PDE.n.nu is not defined. Defaulting to zero');
    psize.n.nu=0;
end

if ~isfield(psize.n,'nw')
    disp('Warning: PDE.n.nw is not defined. Defaulting to zero');
    psize.n.nw=0;
end

if ~isfield(psize.n,'no')
    disp('Warning: PDE.n.no is not defined. Defaulting to zero');
    psize.n.no=0;
end

if ~isfield(psize.n,'noo')
    disp('Warning: PDE.n.noo is not defined. Defaulting to zero');
    psize.n.noo=0;
end

if ~isfield(psize.n,'nro')
    disp('Warning: PDE.n.nro is not defined. Defaulting to zero');
    psize.n.nro=0;
end

if isfield(psize.PDE,'Bpw_nonpol')
    uinput.Bpw_nonpol=psize.PDE.Bpw_nonpol;
end

%  Check disturbance setup

if (psize.n.nw>0)
    if ~isfield(uinput,'w')
        disp('Warning: nw is greater than zero, but user-defiened w inputs are not provided. Defaulting PDE.nw to zero.');
        psize.n.nw=0;
    else
    % Make sure uinput.w is a cell array
    if ~iscell(uinput.w)
        uinput.w=num2cell(uinput.w);
    end
    if (size(uinput.w,2)<psize.n.nw)
        disp('Warning: Size of user-defined disturbances is smaller than PDE.n.nw. Initializing missing inputs to zero');
        uinput.w{1,size(uinput.w,2)+1:psize.n.nw}=0*st;
    end
    end % if ~isfield(uinput,'w')
    end % psize.n.nw>0

    %  Check control inputs setup

if (psize.n.nu>0)
    if ~isfield(uinput,'u')
        disp('Warning: nu is greater than zero, but user-defiened u inputs are not provided. Defaulting PDE.nu to zero.');
        psize.n.nu=0;
    else
    % Make sure uinput.u is a cell array
    if ~iscell(uinput.u)
        uinput.u=num2cell(uinput.u);
    end
    if (size(uinput.u,2)<psize.n.nu)
        disp('Warning: Size of user-defined disturbances is smaller than PDE.n.nw. Initializing missing inputs to zero');
        uinput.u{1,size(uinput.u,2)+1:psize.n.nu}=0*st;
    end
    end % if ~isfield(uinput,'u')
    end % psize.n.nu>0


if ~isfield(psize,'dom')
    disp('PDE domain is not defined. Defaulting to [-1, 1]');
    psize.dom=[-1 1];
end

if (psize.dom(1)==psize.dom(2))
    disp('Warning: left and right ends of the domain are the same. Defaulting domain to [-1, 1]');
    psize.dom=[-1 1];
end

uinput.a=psize.dom(1);
uinput.b=psize.dom(2);

structure=psize;
clear psize;

psize.nu=structure.n.nu;
psize.nw=structure.n.nw;
psize.no=structure.n.no;
psize.nro=structure.n.nro;
psize.noo=structure.n.noo;
psize.N=opts.N;
psize.n=structure.n.n_pde;
psize.dim=1;

ns = sum(psize.n);

 if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.no+psize.n)
            disp('Warning: number of exact solutions provided does not match the number of PDE states');
            disp('Defaulting opts.ifexact to false');
            opts.ifexact=false;
        end % ~isfield(uinput,'exact')
    end % opts.ifexact

    % Check initial conditions

    if ~isfield(uinput,'ic')
    disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
    disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
    uinput.ic.PDE(1:ns)=sin(pi*sx);
    uinput.ic.ODE(1:psize.no)=0;
end

if ~isfield(uinput.ic,'PDE')
    disp('Warning: PDE initial conditions are not defined. Defaulting to sinusoidal');
    uinput.ic.PDE(1:ns)=sin(pi*sx);
elseif ~isempty(symvar(sym(uinput.ic.PDE))) && any(~ismember(symvar(uinput.ic.PDE),{'sx'}))
    error('Initial conditions for PDE must be symbolic expressions in sx');
end

if (psize.no>0)
    if ~isfield(uinput.ic,'ODE')
        disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
        uinput.ic.ODE(1:psize.no)=0;
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

if (psize.no>0)
    if(size(uinput.ic.ODE,2)<psize.no)
        disp('Warning: Number of initial conditions on ODE states is less than the number of ODE states');
        disp('Defalting the rest to zero');
        uinput.ic.ODE(size(uinput.ic.PDE,2)+1:psize.no)=0;
    end

    if(size(uinput.ic.ODE,2)>psize.no)
        disp('Warning: Number of initial conditions on ODE states is greater than the number of ODE states');
        disp('Defaulting all ODE initial conditions to zero');
        clear uinput.ic.ODE;
        uinput.ic=rmfield(uinput.ic,'ODE');
        uinput.ic.ODE(1:psize.no)=0.;
    end
end


end

