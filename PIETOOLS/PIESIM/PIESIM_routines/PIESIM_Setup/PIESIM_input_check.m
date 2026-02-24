%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_input_check.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [structure, uinput, opts, psize]=PIESIM_input_check(varargin)
% Check if necessary inputs are defined to start the simulations
% notE: All other variables will be checked in PIETOOLS converter
% Inputs:
% varargin - variable number of arguments (between 1 and 4)
% Required:
% 1) varargin{1}: data structure of the proglem: PDE, DDE or PIE
% PIE structure of the problem specifies PI operators, T,Tu,Tw, A, Bi, Ci, Dij as fields
% Optional:
% 2) varargin{2}: opts - options for simulation parameters. If empty or incomplete, will be
% set to default values
% 3) varargin{3}: uinput - user-defined boundary inputs, forcing and initial
% conditions. If empty or incomplete, will be set to default values
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
% of PDE. Enhanced the input check for terms format. Terms format has n0w
% its own structure ('PDT'), while batch format has 'PDE' structure.
% Renamed uinput.B21_nonpol to uinput.Bpw_nonpol
% DJ, 10/10/2022 - Update to new terms format. The function can still
% process batch format and old terms format, but these functions have been
% added to separate subroutines.
%   Also set values for no0, nr0 in case of type=='DDE'
% DJ, 12/16/2024: Change default variables to (s1,s1_dum) in ODE case;
% DJ, 12/28/2024: Add spatial domain to "uinput" in 2D case;
% YP, 06/03/2025: Fixed sensing of observed and regulated outputs during coupling with PIE
% YP, 12/29/2025: Updated disturbance and control input check
% YP, 1/27/2026: Added support for infinite-dimensional outputs in 2D
% YP, 2/17/2026: Changed all initial conditions to uinput.ic format

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
    psize.dim = 1;

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
psize.n0 = sum(x_tab(~any(x_tab(:,3:2+PDE.dim),2),2));
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

    % Split vector-valued state components into separate scalar-valued ones
    psize.x_tab = zeros(sum(x_tab(:,2)),size(x_tab,2));
    r_strt = 0;
    for i=1:size(x_tab,1)
        nr = x_tab(i,2);    % Check the size of state component i
        psize.x_tab(r_strt+(1:nr),:) = repmat(x_tab(i,:),nr,1);
        r_strt = r_strt+nr;
    end
    psize.dim=PDE.dim;

%    x_order_PDE = x_order(psize.n0+1:end);
    n0=psize.n0;
    if PDE.dim==1   
    ns = sum(psize.n,'all');
    else
        nx=sum(psize.nx,'all');
        ny=sum(psize.ny,'all');
        n2=sum(psize.n,'all');
    ns = nx+ny+n2;
    end

    % Establish sizes of inputs and outputs
    psize.nw = sum(PDE.w_tab(:,2)); % number of disturbances
    psize.nu = sum(PDE.u_tab(:,2)); % number of control inputs

   % Define sizes for inputs and outputs
    % Determine which inpus and outputs vary in either x, y, or both

     if (PDE.dim==1)
    % Check if the problem size is defined correctly
     if length(opts.N)>1
         disp('Warning: opts.N was defined with more than one entry. The second entry will be ignored.');
         opts.N=opts.N(1);
     end

    psize.nwx = sum(PDE.w_tab(:,3) == 1);
    psize.nwy = 0;
    psize.nw2 = 0;
    psize.nw0 = psize.nw-(psize.nwx+psize.nwy+psize.nw2);

    psize.nux = sum(PDE.u_tab(:,3) == 1);
    psize.nuy = 0;
    psize.nu2 = 0;
    psize.nu0 = psize.nu-(psize.nux+psize.nuy+psize.nu2);

    psize.nrx=sum(PDE.z_tab(:,3)); % number of infinite-dimensional regulated outputs
    psize.nr0 = sum(PDE.z_tab(:,2))-psize.nrx; % number of finite-dimensional regulated outputs
    psize.nox = sum(PDE.y_tab(:,3)); % number of infinite-dimensional observed outputs
    psize.no0 = sum(PDE.y_tab(:,2))-psize.nox; % number of finite-dimensional observed outputs
    % % 
     else
 
    % Check if the problem size is defined correctly
    if (length(opts.N)==1)
        opts.N(2)=opts.N(1);
    elseif size(opts.N,1)>1
        opts.N=opts.N';
    end

    psize.nwx = sum(PDE.w_tab(:,3) == 1 & PDE.w_tab(:,2+PDE.dim) == 0);
    psize.nwy = sum(PDE.w_tab(:,3) == 0 & PDE.w_tab(:,2+PDE.dim) == 1);
    psize.nw2 = sum(PDE.w_tab(:,3) == 1 & PDE.w_tab(:,2+PDE.dim) == 1);
    psize.nw0 = psize.nw-(psize.nwx+psize.nwy+psize.nw2);

    psize.nux = sum(PDE.u_tab(:,3) == 1 & PDE.u_tab(:,2+PDE.dim) == 0);
    psize.nuy = sum(PDE.u_tab(:,3) == 0 & PDE.u_tab(:,2+PDE.dim) == 1);
    psize.nu2 = sum(PDE.u_tab(:,3) == 1 & PDE.u_tab(:,2+PDE.dim) == 1);
    psize.nu0 = psize.nu-(psize.nux+psize.nuy+psize.nu2);

    psize.nrx = sum(PDE.z_tab(:,3) == 1 & PDE.z_tab(:,2+PDE.dim) == 0);
    psize.nry = sum(PDE.z_tab(:,3) == 0 & PDE.z_tab(:,2+PDE.dim) == 1);
    psize.nr2 = sum(PDE.z_tab(:,3) == 1 & PDE.z_tab(:,2+PDE.dim) == 1);
    psize.nr0 = sum(PDE.z_tab(:,2))-(psize.nrx+psize.nry+psize.nr2);

    psize.nox = sum(PDE.y_tab(:,3) == 1 & PDE.y_tab(:,2+PDE.dim) == 0);
    psize.noy = sum(PDE.y_tab(:,3) == 0 & PDE.y_tab(:,2+PDE.dim) == 1);
    psize.no2 = sum(PDE.y_tab(:,3) == 1 & PDE.y_tab(:,2+PDE.dim) == 1);
    psize.no0 = sum(PDE.y_tab(:,2))-(psize.nox+psize.noy+psize.no2);
     end % PDE.dim

    % Define problem size for discretization
    
    psize.N = opts.N;

    if ~isfield(uinput,'ic')
        if n0>0
            disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
        end
        if ns>0
            disp('Warning: PDE initial conditions are not defined. Defaulting to zero');
        end
        uinput.ic(1:n0)=0;
        if (PDE.dim==1)
        uinput.ic(n0+1:n0+ns)=0;
        elseif (PDE.dim==2)
            uinput.ic(n0+1:n0+nx)=0;
            uinput.ic(n0+nx+1:n0+nx+ny)=0;
            uinput.ic(n0+nx+ny+1:n0+n)=0;
        end
   
    else

    ireorder=0;

  
        if ~isempty(symvar(sym(uinput.ic))) && any(~ismember(symvar(uinput.ic),{'sx';'sy'}))
            error('Initial conditions for PDE must be symbolic expressions in sx (and sy for 2D PDEs)');
        elseif length(uinput.ic)>(n0 + ns)
        disp('Warning: Number of initial conditions is greater than the number of states. Extra initial conditions will be ignored.');
        uinput.ic(n0+ns+1:length(uinput.ic))=[];
        elseif length(uinput.ic)<(n0 + ns)
        disp('Warning: Number of initial conditions is less than the number of states. Defaulting the rest to zer.');
        uinput.ic(length(uinput.ic)+1:n0+ns)=sym(0);
        end
        if ~issorted(x_order)
        uinput.ic = uinput.ic(x_order); % Reorder initial conditions to match new ordering of state components.
        ireorder=1;
        disp('Initial conditions have been reordered to correspond to new ordering');
        end
        if (n0>0)
            ic_ode=uinput.ic(1:n0);
            if isa(ic_ode, 'sym')
        if hasSymType(ic_ode,'variable')
            disp('Warning: initial conditions for ODE must be scalar constants. Defaulting to zero');
            uinput.ic(1:n0)=0;
        end
            end
        end
    end % ~isfield(uinput,'ic')

    if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.n0+ns)
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
    psize.n0=PIE.T.dim(1,1);
    nx = PIE.T.dim(2,1);
    psize.N=opts.N;
    psize.n=[0 nx];
    %------------------------
    % Size of finite-dimensional and infinite-dimensional outputs
    psize.nr0 = size(PIE.C1.Q1,1);
    psize.no0 = size(PIE.C2.Q1,1); 
    psize.nrx = size(PIE.C1,1)-psize.nr0;
    psize.nox = size(PIE.C2,1)-psize.no0;  

    psize.nu=size(PIE.Tu,2); % number of disturbances
    psize.nw=size(PIE.Tw,2); % number of control inputs

    % Check if the problem size is defined correctly
    if length(opts.N)>1
         disp('Warning: opts.N was defined with more than one entry. The second entry will be ignored.');
         opts.N=opts.N(1);
     end

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

    n0=psize.n0;
    if ~isfield(uinput,'ic')
        if (n0>0)
        disp('Initial conditions on ODE states are not defined. Defaulting to zero');
        uinput.ic(1:n0)=0;
        end
        if nx>0
        disp('Initial conditions on history are not defined. Defaulting to zero');
        uinput.ic(n0+1:n0+nx)=0;
        end
    elseif(length(uinput.ic)<n0+nx)
         disp('Warning: Number of initial conditions is less than the number of states');
         disp('Defalting the rest to zero');
         uinput.ic(length(uinput.ic)+1:n0+nx)=0;
    elseif(length(uinput.ic)>n0+nx)
         disp('Warning: Number of initial conditions on DDE states is greater than the number of DDE states');
         disp('Extra initial conditions will be ignored');
         uinput.ic(n0+ns+1:length(uinput.ic))=[];
    end

        if (n0>0)
            ic_ode=uinput.ic(1:n0);
            if isa(ic_ode, 'sym')
        if hasSymType(ic_ode,'variable')
            disp('Warning: initial conditions for ODE must be scalar constants. Defaulting to zero');
            uinput.ic(1:n0)=0;
        end
            end
        end


     if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.n0+ns)
            disp('Warning: number of exact solutions provided does not match the number of PDE states');
            disp('Defaulting opts.ifexact to false');
            opts.ifexact=false;
        end % ~isfield(uinput,'exact')
    end % opts.ifexact


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
    psize.n0=PIE.T.dim(1,1);
    nx = PIE.T.dim(2,1);
    psize.N=opts.N;

     % Get the number of PDE states differentiable up to each order
    if isempty(PIE.x_tab) 
    psize.n=nx;
    else
    psize.n = zeros(1,max(PIE.x_tab(:,end))+1);
    for jj=0:max(PIE.x_tab(:,end))
        psize.n(jj+1) = sum(PIE.x_tab(PIE.x_tab(:,end)==jj & PIE.x_tab(:,3),2));
    end
    end

% Find the number of finite-dimensional and infinite-dimensional outputs
% nr0, no0 counts finite-dimensional outputs
% nrx, nox counts infinite-dimensional outputs

    %------------------------
    psize.nr0 = size(PIE.C1.Q1,1);
    psize.no0 = size(PIE.C2.Q1,1); 
    psize.nrx = size(PIE.C1,1)-psize.nr0;
    psize.nox = size(PIE.C2,1)-psize.no0;  

    psize.nu=size(PIE.Tu,2); % number of disturbances
    psize.nw=size(PIE.Tw,2); % number of control inputs

    % Determine which disturbances vary in either x, y, or both

     if (psize.dim==1)

            % Check if the problem size is defined correctly
    if length(opts.N)>1
         disp('Warning: opts.N was defined with more than one entry. The second entry will be ignored.');
         opts.N=opts.N(1);
     end
    
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

    n0=psize.n0;
    
    if ~isfield(uinput,'ic')
        if (n0>0)
        disp('Initial conditions on ODE states are not defined. Defaulting to zero');
        uinput.ic(1:n0)=0;
        end
        if (nx>0)
        disp('Initial conditions on PIE states  are not defined. Defaulting to zero');
        uinput.ic(n0+1:n0+nx)=0;
        end
    elseif(length(uinput.ic)<n0+nx)
         disp('Warning: Number of initial conditions is less than the number of states');
         disp('Defalting the rest to zero');
         uinput.ic(length(uinput.ic)+1:n0+nx)=0;
    elseif(length(uinput.ic)>n0+nx)
         disp('Warning: Number of initial conditions on PIE states is greater than the number of PIE states');
         disp('Extra initial conditions will be ignored');
         uinput.ic(n0+nx+1:length(uinput.ic))=[];
    end

     if (n0>0)
            ic_ode=uinput.ic(1:n0);
            if isa(ic_ode, 'sym')
        if hasSymType(ic_ode,'variable')
            disp('Warning: initial conditions for ODE must be scalar constants. Defaulting to zero');
            uinput.ic(1:n0)=0;
        end
            end
        end


 if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.n0+nx)
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

if length(opts.N)>1
         disp('Warning: opts.N was defined with more than one entry. The second entry will be ignored.');
         opts.N=opts.N(1);
     end

if ~isfield(PDE,'n0')
    disp('Warning: number of PDE states with differentiability of zero were not defined. Defaulting to zero');
    PDE.n0=0;
end
if ~isfield(PDE,'n1')
    disp('Warning: number of PDE states with differentiability of one is not defined. Defaulting to zero');
    PDE.n1=0;
end
if ~isfield(PDE,'n2')
    disp('Warning: number of PDE states with differentiability of two is not defined. Defaulting to zero');
    PDE.n2=0;
end

psize.n=[PDE.n0 PDE.n1 PDE.n2];
ns = PDE.n0+PDE.n1+PDE.n2;

if ~isfield(PDE,'no')
    disp('Warning: number of ODE states is not defined. Defaulting to zero');
    PDE.no=0;
end
if ~isfield(PDE,'nw')
    disp('Warning: number of disturbances is not defined. Defaulting to zero');
    PDE.nw=0;
end
if ~isfield(PDE,'nu')
    disp('Warning: number of control inputs is not defined. Defaulting to zero');
    PDE.nu=0;
end
if ~isfield(PDE,'nb')
    disp('Warning: nb number is not defined. Defaulting to zero');
    PDE.nb=0;
end
if ~isfield(PDE,'nz')
    disp('Warning: number of regulated outputs is not defined. Defaulting to zero');
    PDE.nz=0;
end
if ~isfield(PDE,'ny')
    disp('Warning: number of observed outputs is not defined. Defaulting to zero');
    PDE.ny=0;
end

% Check initial conditions

n0=PDE.no;

if ~isfield(uinput,'ic')
    if (n0>0)
    disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
    uinput.ic(1:n0)=0;
    end
 if (ns>0)
    disp('Warning: PDE initial conditions are not defined. Defaulting to zero');
    uinput.ic(n0+1:n0+ns)=0;
 end
elseif ~isempty(symvar(sym(uinput.ic))) && any(~ismember(symvar(uinput.ic),{'sx'}))
    error('Initial conditions for PDE must be symbolic expressions in sx');
end

if(length(uinput.ic)<n0+ns)
    disp('Warning: Number of initial conditions is less than the number of states');
    disp('Defalting the rest to zero');
    uinput.ic(length(uinput.ic)+1:n0+ns)=0;
end

if(length(uinput.ic)>n0+ns)
    disp('Warning: Number of initial conditions is greater than the number of states');
    disp('Extra initial conditions will be ignored');
    uinput.ic(n0+ns+1:length(uinput.ic))=[];
end

       if (n0>0)
            ic_ode=uinput.ic(1:n0);
            if isa(ic_ode, 'sym')
        if hasSymType(ic_ode,'variable')
            disp('Warning: initial conditions for ODE must be scalar constants. Defaulting to zero');
            uinput.ic(1:n0)=0;
        end
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

structure=PDE;

% Define problem size for discretization

psize.nu=PDE.nu; % number of control inputs
psize.nw=PDE.nw; % number of disturbances
psize.n0=PDE.no; % number of ODE states
psize.nr0=PDE.nz; % number of regulated outputs
psize.no0=PDE.ny; % number of observed outputs
psize.N=opts.N; % order of disceretization in space
psize.dim=1;


 if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.n0+ns)
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

if length(opts.N)>1
         disp('Warning: opts.N was defined with more than one entry. The second entry will be ignored.');
         opts.N=opts.N(1);
end

if ~isfield(psize.n,'nu')
    disp('Warning: number of control inputs is not defined. Defaulting to zero');
    psize.n.nu=0;
end

if ~isfield(psize.n,'nw')
    disp('Warning: number of disturbances is not defined. Defaulting to zero');
    psize.n.nw=0;
end

if ~isfield(psize.n,'nx')
    disp('Warning: number of ODE state variables is not defined. Defaulting to zero');
    psize.n.nx=0;
end

if ~isfield(psize.n,'nz')
    disp('Warning: number of regulated outputs is not defined. Defaulting to zero');
    psize.n.nz=0;
end

if ~isfield(psize.n,'ny')
    disp('Warning: number of observed outputsis not defined. Defaulting to zero');
    psize.n.ny=0;
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
psize.n0=structure.n.nx;
psize.nr0=structure.n.nz;
psize.no0=structure.n.ny;
psize.N=opts.N;
psize.n=structure.n.n_pde;
psize.dim=1;

ns = sum(psize.n);

 if (opts.ifexact)
        if ~isfield(uinput,'exact')
            disp('Warning: exact solution is not provided. Defaulting opts.ifexact to false');
            opts.ifexact=false;
        elseif (size(uinput.exact,2)~=psize.n0+psize.n)
            disp('Warning: number of exact solutions provided does not match the number of PDE states');
            disp('Defaulting opts.ifexact to false');
            opts.ifexact=false;
        end % ~isfield(uinput,'exact')
    end % opts.ifexact

    % Check initial conditions

    n0=psize.n0;

  if ~isfield(uinput,'ic')
    if (n0>0)
    disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
    uinput.ic(1:n0)=0;
    end
 if (ns>0)
    disp('Warning: PDE initial conditions are not defined. Defaulting to zero');
    uinput.ic(n0+1:n0+ns)=0;
 end
elseif ~isempty(symvar(sym(uinput.ic))) && any(~ismember(symvar(uinput.ic),{'sx'}))
    error('Initial conditions for PDE must be symbolic expressions in sx');
end

if(length(uinput.ic)<n0+ns)
    disp('Warning: Number of initial conditions is less than the number of states');
    disp('Defalting the rest to zero');
    uinput.ic(length(uinput.ic)+1:n0+ns)=0;
end

if(length(uinput.ic)>n0+ns)
    disp('Warning: Number of initial conditions is greater than the number of states');
    disp('Extra initial conditions will be ignored');
    uinput.ic(n0+ns+1:length(uinput.ic))=[];
end

        if (n0>0)
            ic_ode=uinput.ic(1:n0);
            if isa(ic_ode, 'sym')
        if hasSymType(ic_ode,'variable')
            disp('Warning: initial conditions for ODE must be scalar constants. Defaulting to zero');
            uinput.ic(1:n0)=0;
        end
            end
        end



end

