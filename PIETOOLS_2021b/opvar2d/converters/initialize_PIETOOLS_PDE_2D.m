%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_PDE_terms_2D.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDE_out=initialize_PIETOOLS_PDE_2D(PDE)
% Initializes and checks the PDE formulation for 2D PDEs in standardized
% format. Undefined signals are set to length 0 with 0 value
% corresponding parameters are set to zero valued terms of appropriate
% dimension

% Note that functionality for 2D PDEs is limited. No inputs or outputs are
% allowed, and the system has to be specified in a very particular manner.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-Defined Inputs for specifying the dynamics:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: All dimensions (nx,nw,nu,nz,ny,nv,nr) are assumed to be 0 by default.
% Specifically, (all optional)
% PDE.n.nx     -  number of ODE states
% PDE.n.nw     -  number of disturbances
% PDE.n.nu     -  number of controlled inputs
% PDE.n.nz     -  number of regulated outputs
% PDE.n.ny     -  number of observed outputs
% PDE.n.nv     -  number of ODE to PDE interconnection signals
% PDE.n.nr     -  number of PDE to ODE interconnection signals
%
% NOTE: Any matrix with a 0-dimension should be ommitted
%
% PDE Terms (required)
% PDE.n.n_pde  -  matrix of dimensions of all PDE states segregated based on differentiability (Required)
% PDE.dom      -  interval of the domain of the spatial variable - s \in [a,b] (Required)

% PDE.PDE.A  -  structure describing the dynamics of the PDE subsystem
% PDE.PDE.Bpv(s) (optional) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nv       - Distributed effect of ODE interconnection in the domain of the PDE

%%%%%% Boundary conditions (required)
% PDE.BC.Ebb     - opvar2d object specifying boundary conditions

% PDE.ODE.Cv    -  coefficients of ODE in v, nv x no size matrix
% PDE.ODE.Dvw    -  coefficients of disturbance in v, nv x nw size matrix
% PDE.ODE.Dvu    -  coefficients of controlled inputs in v, nv x nu size matrix


% PDE to ODE
% PDE.PDE.Crp - structure specifying integral terms in interconnection signal r, size nr x sum((i+1)*n_pde(i),i=0 to N)
% PDE.PDE.Drv - matrix of dimension nr x nv           - feed through term from v to r
% PDE.PDE.Drb - structure specifying Boundary values terms in interconnection signal r, size nr x 2*sum((i)*n_pde(i),i=1 to N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Dynamics (all optional)
% PDE.n.nx     -  number of ODE states
% PDE.ODE.A      - matrix of dimension no x no                             - ODE dynamics
% PDE.ODE.Bw    -  matrix of dimension no x nw                            - Effect of disturbance on ODE state
% PDE.ODE.Bu    -  matrix of dimension no x nu                            - Effect of input on ODE state
% PDE.ODE.Br    -  matrix of dimension no x nr                            - Effect of interconnection signal r on ODE state
%

% Still missing: PDE.PDE.Drv, PDE.PDE.Bpb, PDE.BC.Ebp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS, DJ - 09_29_2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppress warnings if desired
suppress = (evalin('base','exist(''silent_initialize_pde'',''var'')') && evalin('base','silent_initialize_pde'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for domain (mandatory variables)
if ~isfield(PDE,'dom')
    PDE.dom = [0,1;0,1];
    disp('Warning: PDE domain not defined. Defaulting to [0,1;0,1]');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize proper structure if not initialized
if ~isfield(PDE,'ODE')
    PDE.ODE = struct();
end
if ~isfield(PDE,'BC')
    PDE.BC = struct();
end
if ~isfield(PDE,'n')
    PDE.n = struct();
end
if ~isfield(PDE,'PDE')
    PDE.PDE = struct();
    %PDE.n.n_pde = 0;
    %PDE.n.nv = 0;       PDE.n.nr = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we try to identify number of ODE states
% NOTE: Number of PDE states is a mandatory variable and must be explicitly
% defined. 

% First make a list of all variables that are currently defined in PDE
% object
list = [fieldnames(PDE.n); fieldnames(PDE.ODE); fieldnames(PDE.PDE); fieldnames(PDE.BC)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logVal = isdefined(list,{'nx';'Bu';'Bw';'Br';'Cy';'Cz';'Cv'});

if ~isempty(logVal)
    if isfield(PDE.ODE,'A')
        PDE.n.nx = size(PDE.ODE.A,1);
    elseif logVal==2
        PDE.n.nx = size(PDE.ODE.Bu,1);
    elseif logVal==3
        PDE.n.nx = size(PDE.ODE.Bw,1);
    elseif logVal==4
        PDE.n.nx = size(PDE.ODE.Br,1);
    elseif logVal==5
        PDE.n.nx = size(PDE.ODE.Cy,2);
    elseif logVal==6
        PDE.n.nx = size(PDE.ODE.Cz,2);
    elseif logVal==7
        PDE.n.nx = size(PDE.ODE.Cv,2);
    end
else
    if ~suppress
    disp('Warning: Number of ODE state variables is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.n.nx=0;
end

logVal = isdefined(list,{'nu';'Bu';'Dzu';'Dyu';'Dvu'});
if ~isempty(logVal)
    if logVal==2
        PDE.n.nu = size(PDE.ODE.Bu,2);
    elseif logVal==3
        PDE.n.nu = size(PDE.ODE.Dzu,2);
    elseif logVal==4
        PDE.n.nu = size(PDE.ODE.Dyu,2);
    elseif logVal==5
        PDE.n.nu = size(PDE.ODE.Dvu,2);
    end
else
    if ~suppress
    disp('Warning: Number of inputs is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.n.nu=0;
end

logVal = isdefined(list,{'nw';'Bw';'Dzw';'Dyw';'Dvw'});
if ~isempty(logVal)
    if logVal==2
        PDE.n.nw = size(PDE.ODE.Bw,2);
    elseif logVal==3
        PDE.n.nw = size(PDE.ODE.Dzw,2);
    elseif logVal==4
        PDE.n.nw = size(PDE.ODE.Dyw,2);
    elseif logVal==5
        PDE.n.nw = size(PDE.ODE.Dvw,2);
    end
else
    if ~suppress
    disp('Warning: Number of disturbances is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.n.nw=0;
end

logVal = isdefined(list,{'nz';'Cz';'Dzw';'Dzu';'Dzr'});
if ~isempty(logVal)
    if logVal==2        
        PDE.n.nz = size(PDE.ODE.Cz,1);
    elseif logVal==3
        PDE.n.nz = size(PDE.ODE.Dzw,1);
    elseif logVal==4
        PDE.n.nz = size(PDE.ODE.Dzu,1);
    elseif logVal==5
        PDE.n.nz = size(PDE.ODE.Dzr,1);
    end
else
    if ~suppress
    disp('Warning: Number of regulated outputs is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.n.nz=0;
end

logVal = isdefined(list,{'ny';'Cy';'Dyw';'Dyu';'Dyr'});
if ~isempty(logVal)
    if logVal==2
        PDE.n.ny = size(PDE.ODE.Cy,1);
    elseif logVal==3
        PDE.n.ny = size(PDE.ODE.Dyw,1);
    elseif logVal==4
        PDE.n.ny = size(PDE.ODE.Dyu,1);
    elseif logVal==5
        PDE.n.ny = size(PDE.ODE.Dyr,1);
    end
else
    if ~suppress
    disp('Warning: Number of observed outputs is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.n.ny=0;
end

logVal = isdefined(list,{'nv';'Cv';'Dvw';'Dvu'});
if ~isempty(logVal)
    if logVal==2        
        PDE.n.nv = size(PDE.ODE.Cv,1);
    elseif logVal==3
        PDE.n.nv = size(PDE.ODE.Dvw,1);
    elseif logVal==4
        PDE.n.nv = size(PDE.ODE.Dvu,1);
    end
else
    if ~suppress
    disp('Warning: Number of ODE 2 PDE interconnection signals is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.n.nv=0;
end


logVal = isdefined(list,{'nr';'Br';'Dzr';'Dyr';'Crp';'Drb';'Drv'});
if ~isempty(logVal)
    if logVal==2   
        PDE.n.nr = size(PDE.ODE.Br,2);
    elseif logVal==3   
        PDE.n.nr = size(PDE.ODE.Dzr,2);
    elseif logVal==4
        PDE.n.nr = size(PDE.ODE.Dyr,2);
    elseif logVal==5   
        PDE.n.nr = size(PDE.PDE.Crp,1);
    elseif logVal==6
        PDE.n.nr = size(PDE.PDE.Drb,1);
    elseif logVal==7   
        PDE.n.nr = size(PDE.PDE.Drv,1);
    end
else
    if ~suppress
    disp('Warning: Number of PDE 2 ODE interconnection signals is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.n.nr=0;
end


if ~isfield(PDE.n,'n_pde')
    disp('Warning:Number of PDE states n_pde is not defined. Defaulting to zero');
    PDE.n.n_pde=[];
    N = [1,1];
    np =0; np_all_derivatives=0; nBV=[0;0;0]; nBC = [0;0;0];
else
    if any(size(PDE.n.n_pde)==1)
        PDE.n.n_pde = diag(PDE.n.n_pde);
    end
    if any(any(PDE.n.n_pde - diag(diag(PDE.n.n_pde))))
        error('Only PDEs with equal differentiability in both spatial variables are currently supported')
    end
    N = size(PDE.n.n_pde)-1;
    if any(N>=3)
        error('At most 2nd order PDEs are currently supported')
    end
    % Dimensionality of the PDE, should be 2 in this function
    Ndim = length(N);
    % find total PDE states
    np = sum(PDE.n.n_pde(:));
    % find all possible derivatives length
    np_all_derivatives = (N+1+1).*(N+1)/2;
    % Establish the number of boundary values and boundary conditions
    nBC_op = [sum(sum(((0:N(1))'.*(0:N(2))).*PDE.n.n_pde));    % corner values
            sum(sum(((0:N(1))'.*PDE.n.n_pde)));             % x-edge values
            sum(sum(((0:N(2)).*PDE.n.n_pde)))];             % y-edge values
    nBC = sum(nBC_op);
    nBV = [sum(sum((4*(0:N(1))'.*(0:N(2))).*PDE.n.n_pde));    % corner values
            sum(sum((2*(0:N(1))'.*PDE.n.n_pde)));              % x-edge values
            sum(sum((2*(0:N(2)).*PDE.n.n_pde)))];              % y-edge values
end

% Total number of possible derivatives
Ntot = prod(N+1);
% Object used for linear indexing later on
linsz_N = cumprod(N+1);
linsz_N = [1,linsz_N(1:end-1)];

% Build a table of all possible orders of derivatives of the state
Ncell = num2cell(N);
degvals_cell = cellfun(@(m) (0:m)',Ncell,'uni',0);
degvals_grid = degvals_cell;
[degvals_grid{:}] = ndgrid(degvals_cell{:});
deg_table = cell2mat(cellfun(@(m)m(:),degvals_grid,'uni',0));

% BCdegvals_cell = cellfun(@(m) (0:m-1)',Ncell,'uni',0);
% BCdegvals_grid = BCdegvals_cell;
% [BCdegvals_grid{:}] = ndgrid(BCdegvals_cell{:});
% BCdeg_table = cell2mat(cellfun(@(m)m(:),BCdegvals_grid,'uni',0));


% % Also initialize spatial variables
if isfield(PDE,'vars')
    if all(size(PDE.vars)==[2,2]) && ispvar(PDE.vars)
        PDE.vars = PDE.vars;
        ss1 = PDE.vars(1,1);    tt1 = PDE.vars(1,2);
        ss2 = PDE.vars(2,1);    tt2 = PDE.vars(2,2);
    elseif length(PDE.vars)==2 && ispvar(PDE.vars)
        ss1 = PDE.vars(1);      ss2 = PDE.vars(2);
        pvar tt1 tt2
        PDE.vars = [ss1,tt1;ss2,tt2];
    else
        error('Spatial variables must be specified as 2x2 array of pvars')
    end
else
    pvar ss1 ss2 tt1 tt2;
    PDE.vars = [ss1, tt1; ss2, tt2];
end

% Define number of variable explicitly for further reference
nx = PDE.n.nx;  n_pde = PDE.n.n_pde;
nw = PDE.n.nw;  nu = PDE.n.nu;  nr = PDE.n.nr;
nz = PDE.n.nz;  ny = PDE.n.ny;  nv = PDE.n.nv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %
% % % Define parameters related to ODE dynamics
% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if isfield(PDE.ODE,'Axx') && ~isfield(PDE.ODE,'A')
    PDE.ODE.A = PDE.ODE.Axx;
    PDE.ODE = rmfield(PDE.ODE,'Axx');
end
if ~isfield(PDE.ODE,'A')
    if ~suppress
        disp('A is undefined. Defaulting to zero');
    end
    PDE.ODE.A = zeros(nx);
elseif any(size(PDE.ODE.A)~=[nx,nx])
    disp('A has incorrect dimension. Defaulting to zero');
    PDE.ODE.A = zeros(nx);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if isfield(PDE.ODE,'Bxw') && ~isfield(PDE.ODE,'Bw')
    PDE.ODE.Bw = PDE.ODE.Bxw;
    PDE.ODE = rmfield(PDE.ODE,'Bxw');
end
if ~isfield(PDE.ODE,'Bw')
    if ~suppress
    disp('Bw is undefined. Defaulting to zero');
    end
    PDE.ODE.Bw = zeros(nx,nw);
elseif any(size(PDE.ODE.Bw)~=[nx,nw])
    disp('Bw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Bw = zeros(nx,nw);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if isfield(PDE.ODE,'Bxu') && ~isfield(PDE.ODE,'Bu')
    PDE.ODE.Bu = PDE.ODE.Bxu;
    PDE.ODE = rmfield(PDE.ODE,'Bxu');
end
if ~isfield(PDE.ODE,'Bu')
    if ~suppress
    disp('Bu is undefined. Defaulting to zero');
    end
    PDE.ODE.Bu = zeros(nx,nu);
elseif any(size(PDE.ODE.Bu)~=[nx,nu])
    disp('Bu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Bu = zeros(nx,nu);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if isfield(PDE.ODE,'Bxr') && ~isfield(PDE.ODE,'Br')
    PDE.ODE.Br = PDE.ODE.Bxr;
    PDE.ODE = rmfield(PDE.ODE,'Bxr');
end
if ~isfield(PDE.ODE,'Br')
    if ~suppress
    disp('Br is undefined. Defaulting to zero');
    end
    PDE.ODE.Br = zeros(nx,nr);
elseif any(size(PDE.ODE.Br)~=[nx,nr])
    disp('Br has incorrect dimension. Defaulting to zero');
    PDE.ODE.Br = zeros(nx,nr);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize ODE 2 PDE interconnection related parameters
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if isfield(PDE.ODE,'Cvx') && ~isfield(PDE.ODE,'Cv')
    PDE.ODE.Cv = PDE.ODE.Cvx;
    PDE.ODE = rmfield(PDE.ODE,'Cvx');
end
if ~isfield(PDE.ODE,'Cv')
    if ~suppress
    disp('Cv is undefined. Defaulting to zero');
    end
    PDE.ODE.Cv = zeros(nv,nx);
elseif any(size(PDE.ODE.Cv)~=[nv,nx])
    disp('Cv has incorrect dimension. Defaulting to zero');
    PDE.ODE.Cv = zeros(nv,nx);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dvw')
    if ~suppress
    disp('Dvw is undefined. Defaulting to zero');
    end
    PDE.ODE.Dvw = zeros(nv,nw);
elseif any(size(PDE.ODE.Dvw)~=[nv,nw])
    disp('Dvw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dvw = zeros(nv,nw);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dvu')
    if ~suppress
    disp('Dvu is undefined. Defaulting to zero');
    end
    PDE.ODE.Dvu = zeros(nv,nu);
elseif any(size(PDE.ODE.Dvu)~=[nv,nu])
    disp('Dvu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dvu = zeros(nv,nu);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Regulated output related parameters
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if isfield(PDE.ODE,'Czx') && ~isfield(PDE.ODE,'Cz')
    PDE.ODE.Cz = PDE.ODE.Czx;
    PDE.ODE = rmfield(PDE.ODE,'Czx');
end
if ~isfield(PDE.ODE,'Cz')
    if ~suppress
    disp('Cz is undefined. Defaulting to zero');
    end
    PDE.ODE.Cz = zeros(nz,nx);
elseif any(size(PDE.ODE.Cz)~=[nz,nx])
    disp('Cz has incorrect dimension. Defaulting to zero');
    PDE.ODE.Cz = zeros(nz,nx);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dzw')
    if ~suppress
    disp('Dzw is undefined. Defaulting to zero');
    end
    PDE.ODE.Dzw = zeros(nz,nw);
elseif any(size(PDE.ODE.Dzw)~=[nz,nw])
    disp('Dzw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dzw = zeros(nz,nw);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dzu')
    if ~suppress
    disp('Dzu is undefined. Defaulting to zero');
    end
    PDE.ODE.Dzu = zeros(nz,nu);
elseif any(size(PDE.ODE.Dzu)~=[nz,nu])
    disp('Dzu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dzu = zeros(nz,nu);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dzr')
    if ~suppress
    disp('Dzr is undefined. Defaulting to zero');
    end
    PDE.ODE.Dzr = zeros(nz,nr);
elseif any(size(PDE.ODE.Dzr)~=[nz,nr])
    disp('Dzr has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dzr = zeros(nz,nr);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Observed output related parameters
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if isfield(PDE.ODE,'Cyx') && ~isfield(PDE.ODE,'Cy')
    PDE.ODE.Cy = PDE.ODE.Cyx;
    PDE.ODE = rmfield(PDE.ODE,'Cyx');
end
if ~isfield(PDE.ODE,'Cy')
    if ~suppress
    disp('Cy is undefined. Defaulting to zero');
    end
    PDE.ODE.Cy = zeros(ny,nx);
elseif any(size(PDE.ODE.Cy)~=[ny,nx])
    disp('Cy has incorrect dimension. Defaulting to zero');
    PDE.ODE.Cy = zeros(ny,nx);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dyw')
    if ~suppress
    disp('Dyw is undefined. Defaulting to zero');
    end
    PDE.ODE.Dyw = zeros(ny,nw);
elseif any(size(PDE.ODE.Dyw)~=[ny,nw])
    disp('Dyw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dyw = zeros(ny,nw);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dyu')
    if ~suppress
    disp('Dyu is undefined. Defaulting to zero');
    end
    PDE.ODE.Dyu = zeros(ny,nu);
elseif any(size(PDE.ODE.Dyu)~=[ny,nu])
    disp('Dyu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dyu = zeros(ny,nu);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
if ~isfield(PDE.ODE,'Dyr')
    if ~suppress
    disp('Dyr is undefined. Defaulting to zero');
    end
    PDE.ODE.Dyr = zeros(ny,nr);
elseif any(size(PDE.ODE.Dyr)~=[ny,nr])
    disp('Dyr has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dyr = zeros(ny,nr);
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %

% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %
% % % Define parameters related to interior PDE dynamics
% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% % Define contribution of ODE to PDE signal to PDE dynamics
if ~isfield(PDE.PDE,'Bpv')
    if ~suppress
    disp('Bpv is undefined. Defaulting to zero');
    end
    PDE.PDE.Bpv = polynomial(zeros(np,nv));
elseif any(size(PDE.PDE.Bpv)~=[np,nv])
    disp('Bpv has incorrect dimension. Defaulting to zero');
    PDE.PDE.Bpv = polynomial(zeros(np,nv));
end


% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% % Define contribution of boundary states to PDE dynamics
% Initialize a temporary Bpb regardless of user input
clear tmp
sz_Bpb = 2*np_all_derivatives - 2*(N+1);
linsz_Bpb = [1,cumprod(sz_Bpb(1:end-1))];
initBpb = cell(sz_Bpb);

degval = cell(1,Ndim);
Rval = cell(1,Ndim);
for m1 = 1:Ndim
    degval{m1} = zeros(sz_Bpb(m1)/2,1);
    Rval{m1} = zeros(sz_Bpb(m1)/2,1);
    indx1 = 0;
    for m2=0:N(m1)-1
        indx2 = indx1+N(m1)-m2;
        degval{m1}(indx1+1:indx2) = m2;
        Rval{m1}(indx1+1:indx2) = m2+1:N(m1);
        indx1 = indx2;
    end
end

for m=1:numel(initBpb)
    indcs = cell(1,Ndim);
    [indcs{:}] = ind2sub(2*np_all_derivatives - 2*(N+1),m);
    indcs = cell2mat(indcs);
    
    tmp.delta = floor(2*(indcs-1)./sz_Bpb);
    indcs = mod(indcs-1,sz_Bpb/2)+1;
    
    tmp.D = cell2mat(cellfun(@(m) degval{m}(indcs(m)),num2cell(1:Ndim),'uni',0));
    tmp.Rstate = cell2mat(cellfun(@(m) Rval{m}(indcs(m)),num2cell(1:Ndim),'uni',0));
    indx_r = linsz_N*(tmp.Rstate')+1;
    tmp.coeff = cell(Ntot,1);
    
    for l = 1:Ntot
        tmp.Lstate{l} = deg_table(l,:);
        indx_l = linsz_N*(tmp.Lstate{l}')+1;
        tmp.coeff{l} = polynomial(zeros(n_pde(indx_l),n_pde(indx_r)));
    end    
    initBpb{m} = tmp;
end

if ~isfield(PDE.PDE,'Bpb')
    if ~suppress
        disp('Bpb is undefined. Defaulting to zero');
    end
    PDE.PDE.Bpb = initBpb;
else
    dispVal = 0;
    for m=1:numel(PDE.PDE.Bpb)
        tmp = PDE.PDE.Bpb{m};
        %del = tmp.delta.*sz_Bpb/2;
        if ~isfield(tmp,'D')
            tmp.D = zeros(1,Ndim);
        end
        d = tmp.D;  j = tmp.delta;  r = tmp.Rstate;
        arr = cell2mat(cellfun(@(l) sum(N(l)-1:-1:N(l)-d(l)),num2cell(1:Ndim),'uni',0));
        loc = j.*(np_all_derivatives-N-1) + arr + r;
        loc = linsz_Bpb*(loc'-1)+1;
        indx_r = linsz_N*(r')+1;
        if ~isfield(tmp,'coeff')
            error(['No coefficients specified in PDE.PDE.Bpb{',num2str(m),'}']);
        elseif isdouble(tmp.coeff) || isa(tmp.coeff,'polynomial')
            indx_l = linsz_N*(tmp.Lstate')+1;   % linear index of Lstate
            if any(size(tmp.coeff)~=[n_pde(indx_l),n_pde(indx_r)]) % wrong size, dont change initBpb
                dispVal = 1;
            else % correct term, replace in initA
                initBpb{loc}.coeff{indx_l} = initBpb{loc}.coeff{indx_l} + tmp.coeff;
            end
        elseif iscell(tmp.coeff)
            for l=1:length(tmp.Lstate)
                indx_l = linsz_N*(tmp.Lstate{l}')+1;
                if any(size(tmp.coeff{l})~=[n_pde(indx_l),n_pde(indx_r)]) % wrong size, dont change initA
                    dispVal = 1;
                else % correct term, replace in initA
                    initBpb{loc}.coeff{indx_l} = initBpb{loc}.coeff{indx_l}+tmp.coeff{l};
                end
            end
        else
            error('Coefficients for "PDE.PDE.Bpb{i,j}" must be specified as a cell or as an object of type "double" or "polynomial"')
        end
        % finally, sub PDE.PDE.Bpb with initDpb which has updated values, correct
        % size and all terms
        PDE.PDE.Bpb = initBpb;
        if dispVal
            disp('Some terms in Bpb have incorrect dimension or are missing. Defaulting them to zero');
        end
    end
end

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% % Define contribution of interior PDE state to PDE dynamics
% initialize a temporary A regardless of user input
clear tmp
initA = cell(N+1);
for indx = 1:Ntot
    tmp = struct();
    tmp.D = deg_table(indx,:);
    R_array = deg_table(all(deg_table>=tmp.D,2),:);
    for l = 1:Ntot
        tmp.Lstate{l} = deg_table(l,:);
        for r=1:size(R_array,1)
            tmp.Rstate{r} = R_array(r,:);
            indx_l = linsz_N*(tmp.Lstate{l}')+1;
            indx_r = linsz_N*(tmp.Rstate{r}')+1;
            tmp.coeff{l,r} = polynomial(zeros(n_pde(indx_l),n_pde(indx_r)));
        end
    end
    initA{indx} = tmp;
end

% if A is undefined define
if ~isfield(PDE.PDE, 'A')
    if ~suppress
    disp('PDE dynamics structure A is undefined. Defaulting to zero');
    end
    PDE.PDE.A = initA;
else % if A is defined, check if the dimensions are correct and check if terms are missing
    dispVal = 0;
    for m=1:numel(PDE.PDE.A)      
        tmp = PDE.PDE.A{m};
        if ~isfield(tmp,'D')
            tmp.D = zeros(1,Ndim);
        end
        if ~isfield(tmp,'I')
            tmp.I = zeros(1,Ndim);
        elseif ~all(tmp.I==0)
            error('Integral terms are currently not supported')
        end
        loc = tmp.D*linsz_N' + 1; %linear index associated with degree tmp.D
        if ~isfield(tmp,'coeff')
            error(['No coefficients specified in PDE.PDE.A{',num2str(m),'}']);
        elseif isdouble(tmp.coeff) || isa(tmp.coeff,'polynomial')
            indx_l = linsz_N*(tmp.Lstate')+1;   % linear index of Lstate
            indx_r = linsz_N*(tmp.Rstate')+1;   % linear index of Rstate
            linsz_R = cumprod(N-tmp.D+1);
            linsz_R = [1,linsz_R(1:end-1)];
            indx_R = linsz_R*(tmp.Rstate-tmp.D)'+1; % linear index of Rstate in set of at least D differentiable states
            if any(size(tmp.coeff)~=[n_pde(indx_l),n_pde(indx_r)]) % wrong size, dont change initA
                dispVal = 1;
            else % correct term, replace in initA
                initA{loc}.coeff{indx_l,indx_R} = initA{loc}.coeff{indx_l,indx_R} + tmp.coeff;
            end
        elseif iscell(tmp.coeff)
            for l=1:length(tmp.Lstate)
                for r=1:length(tmp.Rstate)
                    indx_l = linsz_N*(tmp.Lstate{l}')+1;
                    indx_r = linsz_N*(tmp.Rstate{r}')+1;
                    linsz_R = cumprod(N-tmp.D+1);
                    linsz_R = [1,linsz_R(1:end-1)];
                    indx_R = linsz_R*(tmp.Rstate{r}-tmp.D)'+1; % linear index of Rstate in set of at least D differentiable states
                    if any(size(tmp.coeff{l,r})~=[n_pde(indx_l),n_pde(indx_r)]) % wrong size, dont change initA
                        dispVal = 1;
                    else % correct term, replace in initA
                        initA{loc}.coeff{indx_l,indx_R} = ...
                            initA{loc}.coeff{indx_l,indx_R}+tmp.coeff{l,r};
                    end
                end
            end
        else
            error('Coefficients for "PDE.PDE.A{i,j}" must be specified as a cell or as an object of type "double" or "polynomial"')
        end
    end
    % finally, sub PDE.PDE.A with initA which has update values, correct
    % size and missing terms
    PDE.PDE.A = initA;
    if dispVal
            disp('Some of the terms in PDE dynamics structure A have incorrect dimension or are missing. Defaulting them to zero');
    end
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %

% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %
% % % Initialize PDE to ODE signal related parameters
% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% % Define contribution of ODE to PDE signal to PDE to ODE signal
% initialize PDE to ODE signal parameters
if isfield(PDE.PDE,'Drv')
    disp('PDE to ODE feedthrough term is currently not supported in 2D; removing field PDE.ODE.Drv.')
    PDE.PDE = rmfield(PDE.PDE,'Drv');
end

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% % Define contribution of boundary PDE states to PDE to ODE signal
% initialize a temporary Drb regardless of user input
clear tmp
sz_Drb = 2*np_all_derivatives - 2*(N+1);
linsz_Drb = [1,cumprod(sz_Drb(1:end-1))];
initDrb = cell(sz_Drb);
% 
% degval = cell(1,Ndim);
% Rval = cell(1,Ndim);
% for m1 = 1:Ndim
%     degval{m1} = zeros(sz_Drb(m1)/2,1);
%     Rval{m1} = zeros(sz_Drb(m1)/2,1);
%     indx1 = 0;
%     for m2=0:N(m1)-1
%         indx2 = indx1+N(m1)-m2;
%         degval{m1}(indx1+1:indx2) = m2;
%         Rval{m1}(indx1+1:indx2) = m2+1:N(m1);
%         indx1 = indx2;
%     end
% end

for m=1:numel(initDrb)
    indcs = cell(1,Ndim);
    [indcs{:}] = ind2sub(2*np_all_derivatives - 2*(N+1),m);
    indcs = cell2mat(indcs);
    
    tmp.delta = floor(2*(indcs-1)./sz_Drb);
    indcs = mod(indcs-1,sz_Drb/2)+1;
    
    tmp.D = cell2mat(cellfun(@(m) degval{m}(indcs(m)),num2cell(1:Ndim),'uni',0));
    tmp.Rstate = cell2mat(cellfun(@(m) Rval{m}(indcs(m)),num2cell(1:Ndim),'uni',0));
    indx_r = linsz_N*(tmp.Rstate')+1;
    tmp.coeff = polynomial(zeros(nr,n_pde(indx_r)));
    
    initDrb{m} = tmp;
end

if ~isfield(PDE.PDE,'Drb')
    if ~suppress
        disp('Drb is undefined. Defaulting to zero');
    end
    PDE.PDE.Drb = initDrb;
else
    dispVal = 0;
    for m=1:numel(PDE.PDE.Drb)
        tmp = PDE.PDE.Drb{m};
        %del = tmp.delta.*sz_Drb/2;
        if ~isfield(tmp,'D')
            tmp.D = zeros(1,Ndim);
        end
        d = tmp.D;  j = tmp.delta;  r = tmp.Rstate;
        arr = cell2mat(cellfun(@(l) sum(N(l)-1:-1:N(l)-d(l)),num2cell(1:Ndim),'uni',0));
        loc = j.*(np_all_derivatives-N-1) + arr + r;
        loc = linsz_Drb*(loc'-1)+1;
        indx_r = linsz_N*(r')+1;
        if any(size(tmp.coeff)~=[nr,n_pde(indx_r)]) 
            % wrong size, dont change initDrb
            dispVal = 1;
        else % correct term, replace in initDrb
            initDrb{loc}.coeff = initDrb{loc}.coeff + tmp.coeff;
        end
        %finally, sub PDE.PDE.Drb with initDpb which has update values, correct
        %size and missing terms
        PDE.PDE.Drb = initDrb;
        if dispVal
            disp('Some terms in Drb have incorrect dimension or are missing. Defaulting them to zero');
        end
    end
end
    
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% % Define contribution of interior PDE states to PDE to ODE signal
% Initialize a temporary Crp regardless of user input
clear tmp
initCrp = cell(N+1);
for indx = 1:Ntot
    tmp = struct();
    tmp.D = deg_table(indx,:);
    R_array = deg_table(all(deg_table>=tmp.D,2),:);
    for r=1:size(R_array,1)
        tmp.Rstate{r} = R_array(r,:);
        indx_r = linsz_N*(tmp.Rstate{r}')+1;
        tmp.coeff{1,r} = polynomial(zeros(nr,n_pde(indx_r)));
    end
    initCrp{indx} = tmp;
end
if ~isfield(PDE.PDE, 'Crp')
    if ~suppress
        disp('PDE output structure Crp is undefined. Defaulting to zero');
    end
    PDE.PDE.Crp = initCrp;
else % if A is defined, check if the dimensions are correct and check if terms are missing
    dispVal = 0;
    for m=1:numel(PDE.PDE.Crp)      
        tmp = PDE.PDE.Crp{m};
        if ~isfield(tmp,'D')
            tmp.D = zeros(1,Ndim);
        end
%         if ~isfield(tmp,'I')
%             tmp.I = zeros(1,Ndim);
%         elseif ~all(tmp.I==0)
%             error('Integral terms are currently not supported')
%         end
        loc = tmp.D*linsz_N' + 1; %linear index associated with degree tmp.D
        if ~isfield(tmp,'coeff')
            error(['No coefficients specified in PDE.PDE.Crp{',num2str(m),'}']);
        elseif isdouble(tmp.coeff) || isa(tmp.coeff,'polynomial')
            indx_r = linsz_N*(tmp.Rstate')+1;   % linear index of Rstate
            linsz_R = cumprod(N-tmp.D+1);
            linsz_R = [1,linsz_R(1:end-1)];
            indx_R = linsz_R*(tmp.Rstate-tmp.D)'+1; % linear index of Rstate in set of at least D differentiable states
            if any(size(tmp.coeff)~=[nr,n_pde(indx_r)]) % wrong size, dont change initA
                dispVal = 1;
            else % correct term, replace in initA
                initCrp{loc}.coeff{1,indx_R} = initCrp{loc}.coeff{1,indx_R} + tmp.coeff;
            end
        elseif iscell(tmp.coeff)
            for r=1:length(tmp.Rstate)
                indx_r = linsz_N*(tmp.Rstate{r}')+1;
                linsz_R = cumprod(N-tmp.D+1);
                linsz_R = [1,linsz_R(1:end-1)];
                indx_R = linsz_R*(tmp.Rstate{r}-tmp.D)'+1; % linear index of Rstate in set of at least D differentiable states
                if any(size(tmp.coeff{1,r})~=[nr,n_pde(indx_r)]) % wrong size, dont change initCrp
                    dispVal = 1;
                else % correct term, replace in initCrp
                    initCrp{loc}.coeff{1,indx_R} = ...
                        initCrp{loc}.coeff{1,indx_R}+tmp.coeff{1,r};
                end
            end
        else
            error('Coefficients for "PDE.PDE.Crp{i,j}" must be specified as a cell or as an object of type "double" or "polynomial"')
        end
    end
    % finally, sub PDE.PDE.A with initA which has update values, correct
    % size and missing terms
    PDE.PDE.Crp = initCrp;
    if dispVal
        disp('Some of the terms in PDE dynamics structure Crp have incorrect dimension or missing. Defaulting them to zero');
    end
end
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %

% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %
% % % Initialize BC related parameters
% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = % % %

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% Define contribution of ODE signals to boundary conditions
if ~isfield(PDE.BC,'Ebv')
    if ~suppress
    disp('Ebv is undefined. Defaulting to zero');
    end
    PDE.BC.Ebv = zeros(nBC,nv);
elseif any(size(PDE.BC.Ebv)~=[nBC,nv])
    disp('Ebv has incorrect dimension. Defaulting to zero');
    PDE.BC.Ebv = zeros(nBC,nv);
end

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% Define contribution of interior PDE state to boundary conditions
if isfield(PDE.BC,'Ebp')
    disp('Integral boundary terms are currently not supported in 2D; removing the field.')
    PDE.BC = rmfield(PDE.BC,'Ebp');
end

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %
% Define contribution of boundary PDE state to boundary conditions
Ebb_dim0 = [sum(sum(((0:2)'*(0:2)).*n_pde)), sum(sum((Ndim^2*(0:2)'*(0:2)).*n_pde))];
Ebb_dimx = [sum(sum((0:2)'.*n_pde)), sum(sum(Ndim*(0:2)'.*n_pde))];
Ebb_dimy = [sum(sum((0:2).*n_pde)), sum(sum(Ndim*(0:2).*n_pde))];
Ebb_dim = [Ebb_dim0;Ebb_dimx;Ebb_dimy;[0,0]];

if ~isfield(PDE.BC,'Ebb')
    error('No boundary conditions are specified')
elseif ~isa(PDE.BC.Ebb,'opvar2d')
    error('Boundary conditions for 2D PDEs must be specified in terms of opvar2d object Ebb')
elseif ~all(all(PDE.BC.Ebb.dim==Ebb_dim))
    error('Operator Ebb is of incorrect dimensions')
end

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % % %

% And we're done! Just clean up the structure a bit
PDE.n = orderfields(PDE.n,{'nx','n_pde','nw','nu','nz','ny','nr','nv'});
PDE.ODE = orderfields(PDE.ODE,{'A','Bw','Bu','Br','Cz','Cy','Cv','Dzw','Dzu','Dzr','Dyw','Dyu','Dyr','Dvw','Dvu'});
PDE.PDE = orderfields(PDE.PDE,{'A','Bpb','Bpv','Crp','Drb'});
PDE.BC = orderfields(PDE.BC,{'Ebb','Ebv'});

PDE_out = orderfields(PDE,{'dom','n','ODE','PDE','BC','vars'});

end 

function logicalVal = isdefined(defined_list, test_list)
%this function tests if the terms in test_list are present in the
%defined_list and returns the first match location
logicalVal = find(ismember(test_list, defined_list));
if ~isempty(logicalVal)
    logicalVal = logicalVal(1);
end
end