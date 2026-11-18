%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_PDE_terms_legacy.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDE_out=initialize_PIETOOLS_PDE_terms_legacy(PDE)
% initializes and checks the PDE formulation using the legacy term-based
% input format.  undefined signals are set to length 0 with 0 value
% corresponding parameters are set to zero valued terms of appropriate
% dimension

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
% PDE.n.n_pde  -  dimensions of all PDE states segregated based on differentiability (Required)
% PDE.dom    -  interval of the domain of the spatial variable - s \in [a,b] (Required)

% PDE.PDE.A  -  structure describing the dynamics of the PDE subsystem
% PDE.PDE.Bpv(s) (optional) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nv       - Distributed effect of ODE interconnection in the domain of the PDE
% PDE.PDE.Bpb(s) (optional) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x 2sum(i*n_pde(i), i=1,to N )       - Distributed effect of boundary values in the domain of the PDE

%%%%%% Boundary conditions (optional)
% PDE.BC.Ebb     - structure specifying boundary conditions
% PDE.BC.Ebp     - structure specifying integral boundary conditions
% PDE.BC.Ebv     -  matrix of dimension nBC x nv        - effect of ODE interconnection on Boundary Conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interconnection signals (all optional)
% ODE to PDE
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
% PDE.ODE.Bxw    -  matrix of dimension no x nw                            - Effect of disturbance on ODE state
% PDE.ODE.Bxu    -  matrix of dimension no x nu                            - Effect of input on ODE state
% PDE.ODE.Bxr    -  matrix of dimension no x nr                            - Effect of interconnection signal r on ODE state
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regulated outputs (all optional)
% PDE.n.nz     -  number of regulated outputs
% PDE.ODE.Cz     -  matrix of dimension nz x no                            - Effect of ODE states on regulated output
% PDE.ODE.Dzw    -  matrix of dimension nz x nw                            - Effect of disturbance on regulated output (avoid for H2-optimal control problems)
% PDE.ODE.Dzu    -  matrix of dimension nz x nu                            - Effect of control input on regulated output (Recommended for realistic controllers)
% PDE.ODE.Dzr    -  matrix of dimension nz x nr                            - Effect of interconnection signal r on regulated output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observed outputs (all optional)
% PDE.n.ny     -  number of observed outputs
% PDE.ODE.Cy     -  matrix of dimension ny x nx                            - Effect of ODE states on observed output
% PDE.ODE.Dyw    -  matrix of dimension ny x nw                            - Effect of disturbance on observed output (e.g. sensor noise)
% PDE.ODE.Dyu    -  matrix of dimension ny x nu                            - Effect of control input on observed output (rare)
% PDE.ODE.Dyr    -  matrix of dimension ny x nr                            - Effect of interconnection signal r on observed outputs

% ~ NOTE: The workspace is cleaned at the end of the file, but feel free to
% exclude variables from this cleaning process if they are needed for
% anything.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP  - 5_29_2021
% SS - 6/1/2021; modified initA, initCrp etc to allow for repeated terms
% DJ - 12/29/2021: Added option to suppress (less important) warnings
% DJ, 12/07/2024: Use new default vars s1 and s1_dum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppress warnings if desired
suppress = (evalin('base','exist(''silent_initialize_pde'',''var'')') && evalin('base','silent_initialize_pde'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for domain (mandatory variables)
if ~isfield(PDE,'dom')
    PDE.dom = [0,1];
    disp('Warning: PDE domain not defined. Defaulting to [0,1]');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize proper structure if not initialized
if ~isfield(PDE,'ODE')
    PDE.ODE = struct();
    %PDE.n.no = 0;   
    %PDE.n.nz = 0;       PDE.n.nw = 0;
    %PDE.n.ny = 0;       PDE.n.nu = 0;
    %PDE.n.nv = 0;       PDE.n.nr = 0;
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
% Now we try to identify dimensions of inputs, outputs, number of ODE
% states
% NOTE: Number of PDE states is a mandatory variable and must be explicitly
% defined. 

% First make a list of all variables that are currently defined in PDE
% object
list = [fieldnames(PDE.n); fieldnames(PDE.ODE); fieldnames(PDE.PDE); fieldnames(PDE.BC)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logVal = isdefined(list,{'nx';'Bxu';'Bxw';'Bxr';'Cy';'Cz';'Cv'});

if ~isempty(logVal)
    if isfield(PDE.ODE,'A')
        PDE.n.nx = size(PDE.ODE.A,1);
    elseif logVal==2
        PDE.n.nx = size(PDE.ODE.Bxu,1);
    elseif logVal==3
        PDE.n.nx = size(PDE.ODE.Bxw,1);
    elseif logVal==4
        PDE.n.nx = size(PDE.ODE.Bxr,1);
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

logVal = isdefined(list,{'nu';'Bxu';'Dzu';'Dyu';'Dvu'});
if ~isempty(logVal)
    if logVal==2
        PDE.n.nu = size(PDE.ODE.Bxu,2);
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

logVal = isdefined(list,{'nw';'Bxw';'Dzw';'Dyw';'Dvw'});
if ~isempty(logVal)
    if logVal==2
        PDE.n.nw = size(PDE.ODE.Bxw,2);
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


logVal = isdefined(list,{'nr';'Bxr';'Dzr';'Dyr';'Crp';'Drb';'Drv'});
if ~isempty(logVal)
    if logVal==2   
        PDE.n.nr = size(PDE.ODE.Bxr,2);
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
    disp('Warning: Number of PDE states n_pde is not defined. Defaulting to zero');
    PDE.n.n_pde=[];
    N=1;
    np =0; np_all_derivatives=0; nBVs=0; nBC = nBVs/2;
else
    N = length(PDE.n.n_pde)-1;
    % find total PDE states
    np = sum(PDE.n.n_pde(1:N+1));
    % find all possible derivatives length
    % % np_all_derivatives = sum((1:N+1).*PDE.n.n_pde);
    np_all_derivatives = sum(1:N+1);
    % find total possible Boundary values
    nBVs=2*sum((0:N).*PDE.n.n_pde); nBC = nBVs/2;
end

% find specified number of BCs
logVal = isdefined(list,{'Ebb';'Ebv';'Ebp'});
if ~isempty(logVal) && ~(length(PDE.n.n_pde)<2)
    if logVal==1
        nBC = size(PDE.BC.Ebb{end}.coeff,1);
    elseif logVal==2
        nBC = size(PDE.BC.Ebv{end}.coeff,1);
    elseif logVal==3
        nBC = size(PDE.BC.Ebp{end}.coeff,1);
    end
elseif ~exist('nBC','var')
    if ~suppress
    disp('Warning: Number of BCs is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    nBC=0;
end


% Assuming dimensions are now known, initialize the remaining parameters to zero based on
% the dimensions
pvar s1 s1_dum;                                                             % DJ, 12/07/2024;

PDE.vars = [s1,s1_dum];

%
% Define number of variable explicitly for further reference
no = PDE.n.nx;  n_pde = PDE.n.n_pde;
nw = PDE.n.nw;  nu = PDE.n.nu;  nr = PDE.n.nr;
nz = PDE.n.nz;  ny = PDE.n.ny;  nv = PDE.n.nv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize ODE dynamics related parameters
if ~isfield(PDE.ODE,'A')
    if ~suppress
    disp('A is undefined. Defaulting to zero');
    end
    PDE.ODE.A = zeros(no);
elseif any(size(PDE.ODE.A)~=[no,no])
    disp('A has incorrect dimension. Defaulting to zero');
    PDE.ODE.A = zeros(no);
end
if ~isfield(PDE.ODE,'Bxw')
    if ~suppress
    disp('Bxw is undefined. Defaulting to zero');
    end
    PDE.ODE.Bxw = zeros(no,nw);
elseif any(size(PDE.ODE.Bxw)~=[no,nw])
    disp('Bxw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Bxw = zeros(no,nw);
end
if ~isfield(PDE.ODE,'Bxu')
    if ~suppress
    disp('Bxu is undefined. Defaulting to zero');
    end
    PDE.ODE.Bxu = zeros(no,nu);
elseif any(size(PDE.ODE.Bxu)~=[no,nu])
    disp('Bxu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Bxu = zeros(no,nu);
end
if ~isfield(PDE.ODE,'Bxr')
    if ~suppress
    disp('Bxr is undefined. Defaulting to zero');
    end
    PDE.ODE.Bxr = zeros(no,nr);
elseif any(size(PDE.ODE.Bxr)~=[no,nr])
    disp('Bxr has incorrect dimension. Defaulting to zero');
    PDE.ODE.Bxr = zeros(no,nr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize ODE 2 PDE interconnection related parameters
if ~isfield(PDE.ODE,'Cv')
    if ~suppress
    disp('Cv is undefined. Defaulting to zero');
    end
    PDE.ODE.Cv = zeros(nv,no);
elseif any(size(PDE.ODE.Cv)~=[nv,no])
    disp('Cv has incorrect dimension. Defaulting to zero');
    PDE.ODE.Cv = zeros(nv,no);
end
if ~isfield(PDE.ODE,'Dvw')
    if ~suppress
    disp('Dvw is undefined. Defaulting to zero');
    end
    PDE.ODE.Dvw = zeros(nv,nw);
elseif any(size(PDE.ODE.Dvw)~=[nv,nw])
    disp('Dvw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dvw = zeros(nv,nw);
end
if ~isfield(PDE.ODE,'Dvu')
    if ~suppress
    disp('Dvu is undefined. Defaulting to zero');
    end
    PDE.ODE.Dvu = zeros(nv,nu);
elseif any(size(PDE.ODE.Dvu)~=[nv,nu])
    disp('Dvu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dvu = zeros(nv,nu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Regulated output related parameters

if ~isfield(PDE.ODE,'Cz')
    if ~suppress
    disp('Cz is undefined. Defaulting to zero');
    end
    PDE.ODE.Cz = zeros(nz,no);
elseif any(size(PDE.ODE.Cz)~=[nz,no])
    disp('Cz has incorrect dimension. Defaulting to zero');
    PDE.ODE.Cz = zeros(nz,no);
end
if ~isfield(PDE.ODE,'Dzw')
    if ~suppress
    disp('Dzw is undefined. Defaulting to zero');
    end
    PDE.ODE.Dzw = zeros(nz,nw);
elseif any(size(PDE.ODE.Dzw)~=[nz,nw])
    disp('Dzw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dzw = zeros(nz,nw);
end
if ~isfield(PDE.ODE,'Dzu')
    if ~suppress
    disp('Dzu is undefined. Defaulting to zero');
    end
    PDE.ODE.Dzu = zeros(nz,nu);
elseif any(size(PDE.ODE.Dzu)~=[nz,nu])
    disp('Dzu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dzu = zeros(nz,nu);
end
if ~isfield(PDE.ODE,'Dzr')
    if ~suppress
    disp('Dzr is undefined. Defaulting to zero');
    end
    PDE.ODE.Dzr = zeros(nz,nr);
elseif any(size(PDE.ODE.Dzr)~=[nz,nr])
    disp('Dzr has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dzr = zeros(nz,nr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Observed output related parameters

if ~isfield(PDE.ODE,'Cy')
    if ~suppress
    disp('Cy is undefined. Defaulting to zero');
    end
    PDE.ODE.Cy = zeros(ny,no);
elseif any(size(PDE.ODE.Cy)~=[ny,no])
    disp('Cy has incorrect dimension. Defaulting to zero');
    PDE.ODE.Cy = zeros(ny,no);
end
if ~isfield(PDE.ODE,'Dyw')
    if ~suppress
    disp('Dyw is undefined. Defaulting to zero');
    end
    PDE.ODE.Dyw = zeros(ny,nw);
elseif any(size(PDE.ODE.Dyw)~=[ny,nw])
    disp('Dyw has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dyw = zeros(ny,nw);
end
if ~isfield(PDE.ODE,'Dyu')
    if ~suppress
    disp('Dyu is undefined. Defaulting to zero');
    end
    PDE.ODE.Dyu = zeros(ny,nu);
elseif any(size(PDE.ODE.Dyu)~=[ny,nu])
    disp('Dyu has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dyu = zeros(ny,nu);
end
if ~isfield(PDE.ODE,'Dyr')
    if ~suppress
    disp('Dyr is undefined. Defaulting to zero');
    end
    PDE.ODE.Dyr = zeros(ny,nr);
elseif any(size(PDE.ODE.Dyr)~=[ny,nr])
    disp('Dyr has incorrect dimension. Defaulting to zero');
    PDE.ODE.Dyr = zeros(ny,nr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize PDE dynamics related parameters

if ~isfield(PDE.PDE,'Bpv')
    if ~suppress
    disp('Bpv is undefined. Defaulting to zero');
    end
    PDE.PDE.Bpv = polynomial(zeros(np,nv));
elseif any(size(PDE.PDE.Bpv)~=[np,nv])
    disp('Bpv has incorrect dimension. Defaulting to zero');
    PDE.PDE.Bpv = polynomial(zeros(np,nv));
end

% initialize a temporary Bpb regardless of user input
% initBpb = repmat({struct('coeff',{},'Lstate',{},'Rstate',{},'D',{},'delta',{})}, 1,(N+1)*(2*np_all_derivatives-2*(N+1)))
initBpb = cell(1,(N+1)*(2*np_all_derivatives-2*(N+1)));
m=1;
for i=0:N
    for j=0:1
        for k=0:N-1   
            for l=k+1:N
                tmp.coeff = polynomial(zeros(n_pde(i+1),n_pde(l+1)));
                tmp.Lstate = i;
                tmp.Rstate = l;
                tmp.D = k;
                tmp.delta = j;
                initBpb{m} = tmp;
                m=m+1;
            end
        end
    end
end
if ~isfield(PDE.PDE,'Bpb')
    if ~suppress
    disp('Bpb is undefined. Defaulting to zero');
    end
    PDE.PDE.Bpb = initBpb;
else
    dispVal = 0;
    for m=1:length(PDE.PDE.Bpb)
        tmp = PDE.PDE.Bpb{m};
        i = tmp.Lstate; j = tmp.delta; l = tmp.Rstate;
        if ~isfield(tmp,'D')
            tmp.D = 0;
        end
        k = tmp.D; 
        loc = i*(2*np_all_derivatives-2*(N+1)) + j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l; 
        if any(size(tmp.coeff)~=[n_pde(tmp.Lstate+1),n_pde(tmp.Rstate+1)]) 
            % wrong size, dont change initBpb
            dispVal = 1;
        else % correct term, replace in initBpb
            initBpb{loc}.coeff = initBpb{loc}.coeff+tmp.coeff;
        end
    end
    %finally, sub PDE.PDE.Bpb with initBpb which has update values, correct
    %size and missing terms
    PDE.PDE.Bpb = initBpb;
    if dispVal
        disp('Some terms in Bpb have incorrect dimension or are missing. Defaulting them to zero');
    end
end


% initialize a temporary A regardless of user input
clear tmp
initA = cell(1,3*(N+1)*np_all_derivatives);
m=1;
for i=0:2
    for j=0:N
        for k=0:N   
            for l=k:N
                tmp.coeff = polynomial(zeros(n_pde(j+1),n_pde(l+1)));
                tmp.Lstate = j;
                tmp.Rstate = l;
                tmp.D = k;
                tmp.I = i;
                initA{m} = tmp;
                m=m+1;
            end
        end
    end
end

% if A is undefined define
if ~isfield(PDE.PDE, 'A')
    if ~suppress
    disp('PDE dynamics structure A is undefined. Defaulting to zero');
    end
    PDE.PDE.A = initA;
else % if A is defined, check if the dimensions are correct and check if terms are missing
    dispVal = 0;
    for m=1:length(PDE.PDE.A)
        tmp = PDE.PDE.A{m};
        j = tmp.Lstate; l = tmp.Rstate;
        if ~isfield(tmp,'D')
            tmp.D = 0;
        end
        if ~isfield(tmp,'I')
            tmp.I = 0;
        end
        k = tmp.D; i = tmp.I;
        loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; % find this
        if any(size(tmp.coeff)~=[n_pde(tmp.Lstate+1),n_pde(tmp.Rstate+1)]) % wrong size, dont change initA
            dispVal = 1;
        else % correct term, replace in initA
            initA{loc}.coeff = initA{loc}.coeff+tmp.coeff;
        end
    end
    %finally, sub PDE.PDE.A with initA which has update values, correct
    %size and missing terms
    PDE.PDE.A = initA;
    if dispVal
            disp('Some of the terms in PDE dynamics structure A have incorrect dimension or are missing. Defaulting them to zero');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize PDE to ODE signal parameters
if ~isfield(PDE.PDE,'Drv')
    if ~suppress
    disp('Drv is undefined. Defaulting to zero');
    end
    PDE.PDE.Drv = zeros(nr,nv);
elseif any(size(PDE.PDE.Drv)~=[nr,nv])
    disp('Drv has incorrect dimension. Defaulting to zero');
    PDE.PDE.Drv = zeros(nr,nv);
end

% initialize a temporary Drb regardless of user input
clear tmp
initDrb = cell(1,2*np_all_derivatives - 2*(N+1));
m=1;
for j=0:1
    for k=0:N-1   
        for l=k+1:N
            tmp.coeff = zeros(nr,n_pde(l+1));
            tmp.Rstate = l;
            tmp.D = k;
            tmp.delta = j;
            initDrb{m} = tmp;
            m=m+1;
        end
    end
end

if ~isfield(PDE.PDE,'Drb')
    if ~suppress
    disp('Drb is undefined. Defaulting to zero');
    end
    PDE.PDE.Drb = initDrb;
else
    dispVal = 0;
    for m=1:length(PDE.PDE.Drb)
        tmp = PDE.PDE.Drb{m};
        j = tmp.delta; l = tmp.Rstate;
        if ~isfield(tmp,'D')
            tmp.D = 0;
        end
        k = tmp.D;
        loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l;
        if any(size(tmp.coeff)~=[nr,n_pde(tmp.Rstate+1)]) 
            % wrong size, dont change initBpb
            dispVal = 1;
        else % correct term, replace in initBpb
            initDrb{loc}.coeff = initDrb{loc}.coeff+tmp.coeff;
        end
    end
    %finally, sub PDE.PDE.Drb with initDpb which has update values, correct
    %size and missing terms
    PDE.PDE.Drb = initDrb;
    if dispVal
        disp('Some terms in Drb have incorrect dimension or are missing. Defaulting them to zero');
    end
end

% initialize a temporary Crp regardless of user input
clear tmp
initCrp = cell(1,np_all_derivatives);
m=1;

for k=0:N   
    for l=k:N
        tmp.coeff = polynomial(zeros(nr,n_pde(l+1)));
        tmp.Rstate = l;
        tmp.D = k;
        initCrp{m} = tmp;
        m=m+1;
    end
end
% if Crp is undefined define
if ~isfield(PDE.PDE, 'Crp')
    if ~suppress
    disp('PDE output structure Crp is undefined. Defaulting to zero');
    end
    PDE.PDE.Crp = initCrp;
else % if Crp is defined, check if the dimensions are correct and check if terms are missing
    dispVal = 0;
    for m=1:length(PDE.PDE.Crp)
        tmp = PDE.PDE.Crp{m};
        l = tmp.Rstate;
        if ~isfield(tmp,'D')
            tmp.D = 0;
        end
        k = tmp.D;
        loc = sum(N:-1:N-k+1) + l + 1; % find this
        if any(size(tmp.coeff)~=[nr,n_pde(tmp.Rstate+1)]) % wrong size, dont change initA
            dispVal = 1;
        else % correct term, replace in initA
            initCrp{loc}.coeff = initCrp{loc}.coeff+tmp.coeff;
        end
    end
    %finally, sub PDE.PDE.Crp with initCrp which has update values, correct
    %size and missing terms
    PDE.PDE.Crp = initCrp;
    if dispVal
        disp('Some of the terms in PDE output structure Crp have incorrect dimension or are missing. Defaulting them to zero');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize BC related parameters

if ~isfield(PDE.BC,'Ebv')
    if ~suppress
    disp('Ebv is undefined. Defaulting to zero');
    end
    PDE.BC.Ebv = zeros(nBC,nv);
elseif any(size(PDE.BC.Ebv)~=[nBC,nv])
    disp('Ebv has incorrect dimension. Defaulting to zero');
    PDE.BC.Ebv = zeros(nBC,nv);
end

% initialize a temporary Ebp regardless of user input
clear tmp
initEbp = cell(1,np_all_derivatives);
m=1;
for k=0:N   
    for l=k:N
        tmp.coeff = polynomial(zeros(nBC,n_pde(l+1)));
        tmp.Rstate = l;
        tmp.D = k;
        initEbp{m} = tmp;
        m=m+1;
    end
end
% if Ebp is undefined define
if ~isfield(PDE.BC, 'Ebp')
    if ~suppress
    disp('PDE BC structure Ebp is undefined. Defaulting to zero');
    end
    PDE.BC.Ebp = initEbp;
else % if Ebp is defined, check if the dimensions are correct and check if terms are missing
    dispVal = 0;
    for m=1:length(PDE.BC.Ebp)
        tmp = PDE.BC.Ebp{m};
        l = tmp.Rstate;
        if ~isfield(tmp,'D')
            tmp.D = 0;
        end
        k = tmp.D;
        loc = sum(N:-1:N-k+1) + l + 1;
        if any(size(tmp.coeff)~=[nBC,n_pde(tmp.Rstate+1)]) % wrong size, dont change initA
            dispVal = 1;
        else % correct term, replace in initEbp
            initEbp{loc}.coeff = initEbp{loc}.coeff+tmp.coeff;
        end
    end
    %finally, sub PDE.BC.Ebp with initCrp which has update values, correct
    %size and missing terms
    PDE.BC.Ebp = initEbp;
    if dispVal
        disp('Some of the terms in PDE BC structure Ebp have incorrect dimension or are missing. Defaulting them to zero');
    end
end



% initialize a temporary Ebb regardless of user input
clear tmp
initEbb = cell(1,2*np_all_derivatives - 2*(N+1));
m=1;

for j=0:1
    for k=0:N-1   
        for l=k+1:N
            tmp.coeff = zeros(nBC,n_pde(l+1));
            tmp.Rstate = l;
            tmp.D = k;
            tmp.delta = j;
            initEbb{m} = tmp;
            m=m+1;
        end
    end
end

if ~isfield(PDE.BC,'Ebb')
    if ~suppress
    disp('Ebb is undefined. Defaulting to zero');
    end
    PDE.BC.Ebb = initEbb;
else
    dispVal = 0;
    for m=1:length(PDE.BC.Ebb)
        tmp = PDE.BC.Ebb{m};
        j = tmp.delta; l = tmp.Rstate;
        if ~isfield(tmp,'D')
            tmp.D = 0;
        end
        k = tmp.D; 
        loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l; 
        if any(size(tmp.coeff)~=[nBC,n_pde(tmp.Rstate+1)]) 
            % wrong size, dont change initBpb
            dispVal = 1;
        else % correct term, replace in initEbb
            initEbb{loc}.coeff = initEbb{loc}.coeff+tmp.coeff;
        end
    end
    %finally, sub PDE.BC.Ebb with initBpb which has update values, correct
    %size and missing terms
    PDE.BC.Ebb = initEbb;
    if dispVal
        disp('Some terms in Ebb have incorrect dimension or are missing. Defaulting them to zero');
    end
end


clear initA initBpb initCrp initDrb initEbp initEbb;
clear no np nr nu nv nw ny nz n_pde nBC nBVs np_all_derivatives
clear i j k l m N
clear list loc logVal dispVal tmp

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