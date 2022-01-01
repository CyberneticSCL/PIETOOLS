%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_PDE_batch.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDE_out=initialize_PIETOOLS_PDE_batch(PDE)
% initializes and checks the PDE formulation using the batch input format. 
% undefined signals are set to length 0 with 0 value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP  - 5_29_2019
%  MMP  - 5_29_2019 - updated to PDE and PIE structures
%  SS - 5_29_2019 - changed all 'exist' conditional statements to 'isfield' or
%  'isdefined' 
% DJ - 12/29/2021: Added option to suppress (less important) warnings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppress warnings if desired
suppress = (evalin('base','exist(''silent_initialize_pde'',''var'')') && evalin('base','silent_initialize_pde'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for domain (mandatory variables)
% if ~exist('a','var')||~exist('b','var') MP 5/21
%     error('Domain not defined');
% end
if ~isfield(PDE,'dom')
    PDE.dom = [0,1];
    disp('Warning: PDE domain not defined. Defaulting to [0,1]');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we try to identify dimensions of inputs, outputs, number of ODE
% states
% NOTE: Number of PDE states is a mandatory variable and must be explicitly
% defined. 

% First make a list of all variables that are currently defined in PDE
% object
list = fieldnames(PDE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logval = isdefined(list,{'nu';'Bu';'B12';'B22';'D12';'D22'});

if ~isempty(logval)
   if logval==2
        PDE.nu = size(PDE.Bu,2);
    elseif logval==3
        PDE.nu = size(PDE.B12,2);
    elseif logval==4
        PDE.nu = size(PDE.B22,2);
    elseif logval==5
        PDE.nu = size(PDE.D12,2);
    elseif logval==6
        PDE.nu = size(PDE.D22,2);
    end
else
    if ~suppress
    disp('Warning: Number of inputs is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.nu=0;
end

logval = isdefined(list,{'nw';'Bw';'B11';'B21';'D11';'D21'});
if ~isempty(logval)
    if logval==2
        PDE.nw = size(PDE.Bw,2);
    elseif logval==3
        PDE.nw = size(PDE.B11,2);
    elseif logval==4
        PDE.nw = size(PDE.B21,2);
    elseif logval==5
        PDE.nw = size(PDE.D11,2);
    elseif logval==6
        PDE.nw = size(PDE.D21,2);
    end
else
    if ~suppress
    disp('Warning: Number of disturbances is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.nw=0;
end

logval = isdefined(list,{'ny';'C2';'Ca2';'Cb2';'Cc2';'D21';'D22';'C20'});
if ~isempty(logval)
    if logval==2
        PDE.ny = size(PDE.C2,1);
    elseif logval==3
        PDE.ny = size(PDE.Ca2,1);
    elseif logval==4
        PDE.ny = size(PDE.Cb2,1);
    elseif logval==5
        PDE.ny = size(PDE.Cc2,1);
    elseif logval==6
        PDE.ny = size(PDE.D21,1);
    elseif logval==7
        PDE.ny = size(PDE.D22,1);
    elseif logval==8
        PDE.ny = size(PDE.C20,1);
    end
else
    if ~suppress
    disp('Warning: Number of observed outputs is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.ny=0;
end

logval = isdefined(list,{'nz';'C1';'Ca1';'Cb1';'Cc1';'D11';'D12';'C10'});
if ~isempty(logval)
    if logval==2
        PDE.nz = size(PDE.C1,1);
    elseif logval==3
        PDE.nz = size(PDE.Ca1,1);
    elseif logval==4
        PDE.nz = size(PDE.Cb1,1);
    elseif logval==5
        PDE.nz = size(PDE.Cc1,1);
    elseif logval==6
        PDE.nz = size(PDE.D11,1);
    elseif logval==7
        PDE.nz = size(PDE.D12,1);
    elseif logval==8
        PDE.nz = size(PDE.C10,1);
    end
else
    if ~suppress
    disp('Warning: Number of regulated outputs is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.nz=0;
end

logval = isdefined(list,{'nx';'A';'E';'Bx';'E0';'Ea';'Eb';'Ec';'B11';'B12'});
if ~isempty(logval)
    if logval==2
        PDE.nx = size(PDE.A,1);
    elseif logval==3
        PDE.nx = size(PDE.E,2);
    elseif logval==4
        PDE.nx = size(PDE.Bx,2);
    elseif logval==5
        PDE.nx = size(PDE.E0,1);
    elseif logval==6
        PDE.nx = size(PDE.Ea,1);
    elseif logval==7
        PDE.nx = size(PDE.Eb,1);
    elseif logval==8
        PDE.nx = size(PDE.Ec,1);
    elseif logval==9
        PDE.nx = size(PDE.B11,1);
    elseif logval==10
        PDE.nx = size(PDE.B12,1);
    end
else
    if ~suppress
    disp('Warning: Number of ODE states is neither defined explicitly nor through parameters. Defaulting to zero');
    end
    PDE.nx=0;
end


if ~isfield(PDE,'n0')
    if ~suppress
    disp('Warning: Number of non-differentiated states n_0 is not defined. Defaulting to zero');
    end
    PDE.n0=0;
end
if ~isfield(PDE,'n1')
    if ~suppress
    disp('Warning: Number of continuously differentiable states n_1 is not defined. Defaulting to zero');
    end
    PDE.n1=0;
end
if ~isfield(PDE,'n2')
    if ~suppress
    disp('Warning: Number of twice continuously differentiable states n_2 is not defined. Defaulting to zero');
    end
    PDE.n2=0;
end

% for convenience unpack dimensions from the PDE object for repeated use
nx = PDE.nx; nw = PDE.nw; nu = PDE.nu; 
n0 = PDE.n0; n1 = PDE.n1; n2 = PDE.n2;
nz = PDE.nz; ny = PDE.ny;
np = PDE.n0+PDE.n1+PDE.n2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if ODE parameters are defined and are of appropriate
% dimensions
if ~isfield(PDE,'E0')
    if ~suppress
    disp('E0 not defined. Defaulting to zero');
    end
    PDE.E0 = zeros(nx, 2*n1+4*n2);
elseif any(size(PDE.E0)~=[nx,2*n1+4*n2])
    disp('E0 has incorrect dimensions. Defaulting to zero');
    PDE.E0 = zeros(nx, 2*n1+4*n2);
end
if ~isfield(PDE,'A')
    if ~suppress
    disp('A is undefined. Defaulting to zero');
    end
    PDE.A = zeros(nx);
elseif any(size(PDE.A)~=[nx,nx])
    disp('A has incorrect dimension. Defaulting to zero');
    PDE.A = zeros(nx);
end
if ~isfield(PDE,'Ea')
    PDE.Ea=zeros(nx,np);
    if ~suppress
    disp('Ea is undefined. Defaulting to zero');
    end
elseif any(size(PDE.Ea)~=[nx,np])
    PDE.Ea=zeros(nx,np);
    disp('Ea has incorrect dimension. Defaulting to zero');
end
if ~isfield(PDE,'Eb')
    PDE.Eb=zeros(nx,n1+n2);
    if ~suppress
    disp('Eb is undefined. Defaulting to zero');
    end
elseif any(size(PDE.Eb)~=[nx,n1+n2])
    PDE.Eb=zeros(nx,n1+n2);
    disp('Eb has incorrect dimension. Defaulting to zero');
end
if ~isfield(PDE,'Ec')
    PDE.Ec=zeros(nx,n2);
    if ~suppress
    disp('Ec is undefined. Defaulting to zero');
    end
elseif any(size(PDE.Ec)~=[nx,n2])
    PDE.Ec=zeros(nx,n2);
    disp('Ec has incorrect dimension. Defaulting to zero');
end
if ~isfield(PDE,'B11')
    PDE.B11=zeros(nx,nw);
    if ~suppress
    disp('B11 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.B11)~=[nx,nw])
    PDE.B11=zeros(nx,nw);
    disp('B11 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'B12')
    PDE.B12=zeros(nx,nu);
    if ~suppress
    disp('B12 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.B12)~=[nx,nu])
    PDE.B12=zeros(nx,nu);
    disp('B12 has incorrect dimension. Defaulting to zero')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if PDE parameters are defined and are of appropriate
% dimensions
if ~isfield(PDE,'A0')
    if ~suppress
    disp('A0 is undefined. Defaulting to zero');
    end
    PDE.A0=zeros(np);
elseif any(size(PDE.A0)~= [np,np])
    disp('A0 has incorrect dimension. Defaulting to zero');
    PDE.A0=zeros(np);
end

if ~isfield(PDE,'A1')
    if ~suppress
    disp('A1 is undefined. Defaulting to zero');
    end
    PDE.A1=zeros(np,n1+n2);
elseif any(size(PDE.A1)~= [np,n1+n2])
    disp('A1 has incorrect dimension. Defaulting to zero');
    PDE.A1=zeros(np,n1+n2);
end

if ~isfield(PDE,'A2')
    if ~suppress
    disp('A2 is undefined. Defaulting to zero');
    end
    PDE.A2=zeros(np,n2);
elseif any(size(PDE.A2)~=[np,n2])
    disp('A2 has incorrect dimension. Defaulting to zero');
    PDE.A2=zeros(np,n2);
end

if ~isfield(PDE,'E')
    if ~suppress
    disp('E is undefined. Defaulting to zero');
    end
    PDE.E = zeros(np,nx);
elseif any(size(PDE.E)~=[np,nx])
    disp('E has incorrect dimension. Defaulting to zero');
    PDE.E = zeros(np,nx);
end

if ~isfield(PDE,'B21')
    PDE.B21=zeros(np,nw);
    if ~suppress
    disp('B21 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.B21)~=[np,nw])
    PDE.B21=zeros(np,nw);
    disp('B21 has incorrect dimension. Defaulting to zero')
end

if ~isfield(PDE,'B22')
    PDE.B22=zeros(np,nu);
    if ~suppress
    disp('B22 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.B22)~=[np,nu])
    PDE.B22=zeros(np,nu);
    disp('B22 has incorrect dimension. Defaulting to zero')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if BCs are defined and are of appropriate
% dimensions. 
% NOTE: B is a mandatory input
if ~isfield(PDE,'B')
    error('B is undefined. Boundary conditions are necessary for PDEs');
elseif any(size(PDE.B)~=[n1+2*n2,2*n1+4*n2])
    error('B has incorrect dimension. Not enough boundary conditions');
end

if ~isfield(PDE,'Bx')
    if ~suppress
    disp('Bx is undefined. Defaulting to zero');
    end
    PDE.Bx = zeros(n1+2*n2,nx);
elseif any(size(PDE.Bx)~=[n1+2*n2,nx])
    disp('Bx has incorrect dimension. Defaulting to zero');
    PDE.Bx = zeros(n1+2*n2,nx);
end

if ~isfield(PDE,'Bw')
    if ~suppress
    disp('Bw is undefined. Defaulting to zero');
    end
    PDE.Bw = zeros(n1+2*n2,nw);
elseif any(size(PDE.Bw)~=[n1+2*n2,nw])
    disp('Bw has incorrect dimension. Defaulting to zero');
    PDE.Bw = zeros(n1+2*n2,nw);
end

if ~isfield(PDE,'Bu')
    if ~suppress
    disp('Bu is undefined. Defaulting to zero');
    end
    PDE.Bu = zeros(n1+2*n2,nu);
elseif any(size(PDE.Bu)~=[n1+2*n2,nu])
    disp('Bu has incorrect dimension. Defaulting to zero');
    PDE.Bu = zeros(n1+2*n2,nu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if regulated output parameters are defined and are of appropriate
% dimensions
if ~isfield(PDE,'C1')
    PDE.C1=zeros(nz,nx);
    if ~suppress
    disp('C1 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.C1)~=[nz,nx])
    PDE.C1=zeros(nz,nx);
    disp('C1 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'C10')
    PDE.C10=zeros(nz,2*n1+4*n2);
    if ~suppress
    disp('C10 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.C10)~=[nz,2*n1+4*n2])
    PDE.C10=zeros(nz,2*n1+4*n2);
    disp('C10 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'Ca1')
    PDE.Ca1=zeros(nz,np);
    if ~suppress
    disp('Ca1 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.Ca1)~=[nz,np])
    PDE.Ca1=zeros(nz,np);
    disp('Ca1 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'Cb1')
    PDE.Cb1=zeros(nz,n1+n2);
    if ~suppress
    disp('Cb1 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.Cb1)~=[nz,n1+n2])
    PDE.Cb1=zeros(nz,n1+n2);
    disp('Cb1 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'Cc1')
    PDE.Cc1=zeros(nz,n2);
    if ~suppress
    disp('Cc1 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.Cc1)~=[nz,n2])
    PDE.Cc1=zeros(nz,n2);
    disp('Cc1 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'D11')
    PDE.D11=zeros(nz,nw);
    if ~suppress
    disp('D11 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.D11)~=[nz,nw])
    PDE.D11=zeros(nz,nw);
    disp('D11 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'D12')
    PDE.D12=zeros(nz,nu);
    if ~suppress
    disp('D12 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.D12)~=[nz,nu])
    PDE.D12=zeros(nz,nu);
    disp('D12 has incorrect dimension. Defaulting to zero')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if sensed output parameters are defined and are of appropriate
% dimensions
if ~isfield(PDE,'C2')
    PDE.C2=zeros(ny,nx);
    if ~suppress
    disp('C2 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.C2)~=[ny,nx])
    PDE.C2=zeros(ny,nx);
    disp('C2 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'C20')
    PDE.C20=zeros(ny,2*n1+4*n2);
    if ~suppress
    disp('C20 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.C20)~=[ny,2*n1+4*n2])
    PDE.C20=zeros(ny,2*n1+4*n2);
    disp('C20 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'Ca2')
    PDE.Ca2=zeros(ny,np);
    if ~suppress
    disp('Ca2 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.Ca2)~=[ny,np])
    PDE.Ca2=zeros(ny,np);
    disp('Ca2 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'Cb2')
    PDE.Cb2=zeros(ny,n1+n2);
    if ~suppress
    disp('Cb2 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.Cb2)~=[ny,n1+n2])
    PDE.Cb2=zeros(ny,n1+n2);
    disp('Cb2 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'Cc2')
    PDE.Cc2=zeros(ny,n2);
    if ~suppress
    disp('Cc2 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.Cc2)~=[ny,n2])
    PDE.Cc2=zeros(ny,n2);
    disp('Cc2 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'D21')
    PDE.D21=zeros(ny,nw);
    if ~suppress
    disp('D21 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.D21)~=[ny,nw])
    PDE.D21=zeros(ny,nw);
    disp('D21 has incorrect dimension. Defaulting to zero')
end
if ~isfield(PDE,'D22')
    PDE.D22=zeros(ny,nu);
    if ~suppress
    disp('D22 is undefined. Defaulting to zero')
    end
elseif any(size(PDE.D22)~=[ny,nu])
    PDE.D22=zeros(ny,nu);
    disp('D22 has incorrect dimension. Defaulting to zero')
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