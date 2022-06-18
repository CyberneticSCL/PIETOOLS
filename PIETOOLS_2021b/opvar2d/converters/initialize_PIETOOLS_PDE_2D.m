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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interconnection signals (all optional)
% ODE to PDE
% PDE.ODE.Cv    -  coefficients of ODE in v, nv x no size matrix

% PDE to ODE
% PDE.PDE.Drb - structure specifying corner boundary values terms in interconnection signal r, size nr x 2*sum((i+j)*n_pde(i,j),i,j=1 to N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Dynamics (all optional)
% PDE.n.nx       -  number of ODE states
% PDE.ODE.A      - matrix of dimension no x no                             - ODE dynamics
% PDE.ODE.Bxr    -  matrix of dimension no x nr                            - Effect of interconnection signal r on ODE state
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS, DJ - 09_29_2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

logVal = isdefined(list,{'nx';'Bxr';'Cv'});
if ~isempty(logVal)
    if isfield(PDE.ODE,'A')
        PDE.n.nx = size(PDE.ODE.A,1);
    elseif logVal==2
        PDE.n.nx = size(PDE.ODE.Bxr,1);
    elseif logVal==3
        PDE.n.nx = size(PDE.ODE.Cv,2);
    end
else
    %disp('Warning:Number of ODE state variables is neither defined explicitly nor through parameters. Defaulting to zero');
    PDE.n.nx=0;
end

logVal = isdefined(list,{'nv';'Cv'});
if ~isempty(logVal)
    if logVal==2        
        PDE.n.nv = size(PDE.ODE.Cv,1);
    end
else
    %disp('Warning:Number of ODE 2 PDE interconnection signals is neither defined explicitly nor through parameters. Defaulting to zero');
    PDE.n.nv=0;
end


logVal = isdefined(list,{'nr';'Bxr';'Crp';'Drv'});
if ~isempty(logVal)
    if logVal==2   
        PDE.n.nr = size(PDE.ODE.Bxr,2);
    elseif logVal==3   
        PDE.n.nr = size(PDE.PDE.Crp,1);
    elseif logVal==4   
        PDE.n.nr = size(PDE.PDE.Drv,1);
    end
else
    %disp('Warning:Number of PDE 2 ODE interconnection signals is neither defined explicitly nor through parameters. Defaulting to zero');
    PDE.n.nr=0;
end


if ~isfield(PDE.n,'n_pde')
    disp('Warning:Number of PDE states n_pde is not defined. Defaulting to zero');
    PDE.n.n_pde=[];
    N = [1,1];
    np =0; np_all_derivatives=0; nBVs=0; nBC = nBVs/2;
else
    if any(size(PDE.n.n_pde)==1)
        PDE.n.n_pde = diag(PDE.n.n_pde);
    end
    if any(PDE.n.n_pde - diag(diag(PDE.n.n_pde)),'all')
        error('Only PDEs with equal x and y differentiability are currently supported')
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


% Assuming dimensions are now known, initialize the remaining parameters to zero based on
% the dimensions
pvar ss1 ss2 tt1 tt2;

PDE.vars = [ss1, tt1; ss2, tt2];

% Define number of variable explicitly for further reference
nx = PDE.n.nx;  n_pde = PDE.n.n_pde;
nr = PDE.n.nr;  nv = PDE.n.nv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize ODE dynamics related parameters
if isfield(PDE.ODE,'Axx') && ~isfield(PDE.ODE,'A')
    PDE.ODE.A = PDE.ODE.Axx;
    PDE.ODE = rmfield(PDE.ODE,'Axx');
end
if ~isfield(PDE.ODE,'A')
    %disp('A is undefined. Defaulting to zero');
    PDE.ODE.A = zeros(nx);
elseif any(size(PDE.ODE.A)~=[nx,nx])
    disp('A has incorrect dimension. Defaulting to zero');
    PDE.ODE.A = zeros(nx);
end
if isfield(PDE.ODE,'Bxw')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Bxw.')
    PDE.ODE = rmfield(PDE.ODE,'Bxw');
end
if isfield(PDE.ODE,'Bxu')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Bxu.')
    PDE.ODE = rmfield(PDE.ODE,'Bxu');
end
if ~isfield(PDE.ODE,'Bxr')
    %disp('Bxr is undefined. Defaulting to zero');
    PDE.ODE.Bxr = zeros(nx,nr);
elseif any(size(PDE.ODE.Bxr)~=[nx,nr])
    disp('Bxr has incorrect dimension. Defaulting to zero');
    PDE.ODE.Bxr = zeros(nx,nr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize ODE 2 PDE interconnection related parameters
if isfield(PDE.ODE,'Cvx') && ~isfield(PDE.ODE,'Cv')
    PDE.ODE.Cv = PDE.ODE.Cvx;
    PDE.ODE = rmfield(PDE.ODE,'Cvx');
end
if ~isfield(PDE.ODE,'Cv')
    %disp('Cv is undefined. Defaulting to zero');
    PDE.ODE.Cv = zeros(nv,nx);
elseif any(size(PDE.ODE.Cv)~=[nv,nx])
    disp('Cv has incorrect dimension. Defaulting to zero');
    PDE.ODE.Cv = zeros(nv,nx);
end
if isfield(PDE.ODE,'Dvw')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dvw.')
    PDE.ODE = rmfield(PDE.ODE,'Dvw');
end
if isfield(PDE.ODE,'Dvu')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dvu.')
    PDE.ODE = rmfield(PDE.ODE,'Dvu');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Regulated output related parameters

if isfield(PDE.ODE,'Cz')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Cz.')
    PDE.ODE = rmfield(PDE.ODE,'Cz');
end
if isfield(PDE.ODE,'Dzw')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dzw.')
    PDE.ODE = rmfield(PDE.ODE,'Dzw');
end
if isfield(PDE.ODE,'Dzu')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dzu.')
    PDE.ODE = rmfield(PDE.ODE,'Dzu');
end
if isfield(PDE.ODE,'Dzr')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dzr.')
    PDE.ODE = rmfield(PDE.ODE,'Dzr');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Observed output related parameters

if isfield(PDE.ODE,'Cy')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Cy.')
    PDE.ODE = rmfield(PDE.ODE,'Cy');
end
if isfield(PDE.ODE,'Dyw')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dyw.')
    PDE.ODE = rmfield(PDE.ODE,'Dyw');
end
if isfield(PDE.ODE,'Dyu')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dyu.')
    PDE.ODE = rmfield(PDE.ODE,'Dyu');
end
if isfield(PDE.ODE,'Dyr')
    disp('Outputs and inputs are currently not supported in 2D; removing the field ODE.Dyr.')
    PDE.ODE = rmfield(PDE.ODE,'Dyr');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize PDE dynamics related parameters

if ~isfield(PDE.PDE,'Bpv')
    %disp('Bpv is undefined. Defaulting to zero');
    PDE.PDE.Bpv = polynomial(zeros(np,nv));
elseif any(size(PDE.PDE.Bpv)~=[np,nv])
    disp('Bpv has incorrect dimension. Defaulting to zero');
    PDE.PDE.Bpv = polynomial(zeros(np,nv));
end


if isfield(PDE.PDE,'Bpb')
    disp('Boundary state contribution to PDE is currently not supported in 2D; removing the field.')
    PDE.PDE = rmfield(PDE.PDE,'Bpb');
end

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
            tmp.coeff{l,r} = zeros(n_pde(indx_l),n_pde(indx_r));
        end
    end
    initA{indx} = tmp;
end

% if A is undefined define
if ~isfield(PDE.PDE, 'A')
    disp('PDE dynamics structure A is undefined. Defaulting to zero');
    PDE.PDE.A = initA;
else % if A is defined, check if the dimensions are correct and check if terms are missing
    dispVal = 0;
    for m=1:length(PDE.PDE.A)      
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
        elseif isdouble(tmp.coeff)
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
                    if any(size(tmp.coeff{l,r})~=[n_pde(indx_l),n_pde(indx_r)]) % wrong size, dont change initA
                        dispVal = 1;
                    else % correct term, replace in initA
                        initA{loc}.coeff{n_pde(indx_l),n_pde(indx_r)} = ...
                            initA{loc}.coeff{n_pde(indx_l),n_pde(indx_r)}+tmp.coeff{l,r};
                    end
                end
            end
        end
    end
    % finally, sub PDE.PDE.A with initA which has update values, correct
    % size and missing terms
    PDE.PDE.A = initA;
    if dispVal
            disp('Some of the terms in PDE dynamics structure A have incorrect dimension or missing. Defaulting them to zero');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize PDE to ODE signal parameters
if isfield(PDE.PDE,'Drv')
    disp('PDE to ODE feedthrough term is currently not supported in 2D; removing field PDE.ODE.Drv.')
    PDE.PDE = rmfield(PDE.PDE,'Drv');
end

% initialize a temporary Drb regardless of user input
clear tmp
sz_Drb = 2*np_all_derivatives - 2*(N+1);
linsz_Drb = [1,cumprod(sz_Drb(1:end-1))];
initDrb = cell(sz_Drb);

degval = cell(1,Ndim);
Rval = cell(1,Ndim);
for m1 = 1:Ndim
    degval{m1} = zeros(sz_Drb(m1)/2,1);
    Rval{m1} = zeros(sz_Drb(m1)/2,1);
    indx1 = 0;
    for m2=0:N(m1)-1
        indx2 = indx1+N(m1)-m2;
        degval{m1}(indx1+1:indx2) = m2;
        Rval{m1}(indx1+1:indx2) = m2+1:N(m1);
        indx1 = indx2;
    end
end

for m=1:numel(initDrb)
    indcs = cell(1,Ndim);
    [indcs{:}] = ind2sub(2*np_all_derivatives - 2*(N+1),m);
    indcs = cell2mat(indcs);
    
    tmp.delta = floor(2*(indcs-1)./sz_Drb);
    indcs = mod(indcs-1,sz_Drb/2)+1;
    
    tmp.D = cell2mat(cellfun(@(m) degval{m}(indcs(m)),num2cell(1:Ndim),'uni',0));
    tmp.Rstate = cell2mat(cellfun(@(m) Rval{m}(indcs(m)),num2cell(1:Ndim),'uni',0));
    indx_r = linsz_N*(tmp.Rstate')+1;
    tmp.coeff = zeros(nr,n_pde(indx_r));
    
    initDrb{m} = tmp;
end

if ~isfield(PDE.PDE,'Drb')
    %disp('Drb is undefined. Defaulting to zero');
    PDE.PDE.Drb = initDrb;
else
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
    

if isfield(PDE.PDE,'Crp')
    disp('Direct PDE state contribution to ODE is currently not supported in 2D; removing the field.')
    PDE.PDE = rmfield(PDE.BC,'Crp');
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize BC related parameters

if isfield(PDE.BC,'Ebv')
    disp('ODE state contribution to BCs is currently not supported in 2D; removing the field.')
    PDE.BC = rmfield(PDE.BC,'Ebv');
end

if isfield(PDE.BC,'Ebp')
    disp('Integral boundary terms are currently not supported in 2D; removing the field.')
    PDE.BC = rmfield(PDE.BC,'Ebp');
end

Ebb_dim0 = [sum(((0:2)'*(0:2)).*n_pde,'all'), sum((Ndim^2*(0:2)'*(0:2)).*n_pde,'all')];
Ebb_dimx = [sum((0:2)'.*n_pde,'all'), sum(Ndim*(0:2)'.*n_pde,'all')];
Ebb_dimy = [sum((0:2).*n_pde,'all'), sum(Ndim*(0:2).*n_pde,'all')];
Ebb_dim = [Ebb_dim0;Ebb_dimx;Ebb_dimy;[0,0]];

if ~isfield(PDE.BC,'Ebb')
    error('No boundary conditions are specified')
elseif ~isa(PDE.BC.Ebb,'opvar2d')
    error('Boundary conditions for 2D PDEs must be specified in terms of opvar2d object Ebb')
elseif ~all(PDE.BC.Ebb.dim==Ebb_dim,'all')
    error('Operator Ebb is of incorrect dimensions')
end


PDE_out=PDE;
end 

function logicalVal = isdefined(defined_list, test_list)
%this function tests if the terms in test_list are present in the
%defined_list and returns the first match location
logicalVal = find(ismember(test_list, defined_list));
if ~isempty(logicalVal)
    logicalVal = logicalVal(1);
end
end