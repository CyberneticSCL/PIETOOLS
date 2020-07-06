%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_PDE.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initializes and checks the PDE formulation 
% undefined signals are set to length 0 with 0 value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for domain (mandatory variables)
if ~exist('a','var')||~exist('b','var')
    error('Domain not defined');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we try to identify dimensions of inputs, outputs, number of ODE
% states
% NOTE: Number of PDE states is a mandatory variable and must be explicitly
% defined. 

if exist('Bu','var')||exist('B12','var')||exist('B22','var')||exist('D12','var')||exist('D22','var')||exist('nu','var')
    if exist('nu','var')
        nu=nu;
    elseif exist('Bu','var')
        nu = size(Bu,2);
    elseif exist('B12','var')
        nu = size(B12,2);
    elseif exist('B22','var')
        nu = size(B22,2);
    elseif exist('D12','var')
        nu = size(D12,2);
    elseif exist('D22','var')
        nu = size(D22,2);
    end
else
    disp('Warning:Number of inputs is neither defined explicitly nor through parameters. Defaulting to zero');
    nu=0;
end
if exist('Bw','var')||exist('B11','var')||exist('B21','var')||exist('D11','var')||exist('D21','var')||exist('nw','var')
    if exist('nw','var')
        nw=nw;
    elseif exist('Bw','var')
        nw = size(Bw,2);
    elseif exist('B11','var')
        nw = size(B11,2);
    elseif exist('B21','var')
        nw = size(B21,2);
    elseif exist('D11','var')
        nw = size(D11,2);
    elseif exist('D21','var')
        nw = size(D21,2);
    end
else
    disp('Warning:Number of disturbances is neither defined explicitly nor through parameters. Defaulting to zero');
    nw=0;
end
if exist('ny','var')||exist('C2','var')||exist('Ca2','var')||exist('Cb2','var')||exist('Cc2','var')||exist('D21','var')||exist('D22','var')
    if exist('ny','var')
        ny=ny;
    elseif exist('C2','var')
        ny = size(C2,1);
    elseif exist('Ca2','var')
        ny = size(Ca2,1);
    elseif exist('Cb2','var')
        ny = size(Cb2,1);
    elseif exist('Cc2','var')
        ny = size(Cc2,1);
    elseif exist('D21','var')
        ny = size(D21,1);
    elseif exist('D22','var')
        ny = size(D22,1);
    end
else
    disp('Warning:Number of observed outputs is neither defined explicitly nor through parameters. Defaulting to zero');
    ny=0;
end
if exist('nz','var')||exist('C1','var')||exist('Ca1','var')||exist('Cb1','var')||exist('Cc1','var')||exist('D11','var')||exist('D12','var')
    if exist('nz','var')
        nz=nz;
    elseif exist('C1','var')
        nz = size(C1,1);
    elseif exist('Ca1','var')
        nz = size(Ca1,1);
    elseif exist('Cb1','var')
        nz = size(Cb1,1);
    elseif exist('Cc1','var')
        nz = size(Cc1,1);
    elseif exist('D11','var')
        nz = size(D11,1);
    elseif exist('D12','var')
        nz = size(D12,1);
    end
else
    disp('Warning:Number of regulated outputs is neither defined explicitly nor through parameters. Defaulting to zero');
    nz=0;
end

if exist('nx','var')||exist('A','var')||exist('E0','var')||exist('Ea','var')||exist('Eb','var')||exist('Ec','var')||exist('B11','var')||exist('B12','var')||exist('Bx','var')||exist('E','var')
    if exist('nx','var')
        nx=nx;
    elseif exist('Bx','var')
        nx = size(Bx,2);
    elseif exist('E','var')
        nx = size(E,2);
    elseif exist('A','var')
        nx = size(A,1);
    elseif exist('E0','var')
        nx = size(E0,1);
    elseif exist('Ea','var')
        nx = size(Ea,1);
    elseif exist('Eb','var')
        nx = size(Eb,1);
    elseif exist('Ec','var')
        nx = size(Ec,1);
    elseif exist('B11','var')
        nx = size(B11,1);
    elseif exist('B12','var')
        nx = size(B12,1);
    end
else
    disp('Warning:Number of ODE states is neither defined explicitly nor through parameters. Defaulting to zero');
    nx=0;
end
if ~exist('n0','var')
    disp('Warning:Number of non-differentiated states n_0 is not defined. Defaulting to zero');
    n0=0;
end
if ~exist('n1','var')
    disp('Warning:Number of continuously differentiable states n_1 is not defined. Defaulting to zero');
    n1=0;
end
if ~exist('n2','var')
    disp('Warning:Number of twice continuously differentiable states n_2 is not defined. Defaulting to zero');
    n2=0;
end
np = n0+n1+n2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if ODE parameters are defined and are of appropriate
% dimensions
if ~exist('E0','var')
    disp('E0 not defined. Defaulting to zero');
    E0 = zeros(nx, 2*n1+4*n2);
elseif any(size(E0)~=[nx,2*n1+4*n2])
    disp('E0 has incorrect dimensions. Defaulting to zero');
    E0 = zeros(nx, 2*n1+4*n2);
end
if ~exist('A','var')
    disp('A is undefined. Defaulting to zero');
    A = zeros(nx);
elseif any(size(A)~=[nx,nx])
    disp('A has incorrect dimension. Defaulting to zero');
    A = zeros(nx);
end
if ~exist('Ea','var')
    Ea=zeros(nx,np);
    disp('Eais undefined. Defaulting to zero');
elseif any(size(Ea)~=[nx,np])
    Ea=zeros(nx,np);
    disp('Ea has incorrect dimension. Defaulting to zero');
end
if ~exist('Eb','var')
    Eb=zeros(nx,n1+n2);
    disp('Eb is undefined. Defaulting to zero');
elseif any(size(Eb)~=[nx,n1+n2])
    Eb=zeros(nx,n1+n2);
    disp('Eb has incorrect dimension. Defaulting to zero');
end
if ~exist('Ec','var')
    Ec=zeros(nx,n2);
    disp('Ec is undefined. Defaulting to zero');
elseif any(size(Ec)~=[nx,n2])
    Ec=zeros(nx,n2);
    disp('Ec has incorrect dimension. Defaulting to zero');
end
if ~exist('B11','var')
    B11=zeros(nx,nw);
    disp('B11 is undefined. Defaulting to zero')
elseif any(size(B11)~=[nx,nw])
    B11=zeros(nx,nw);
    disp('B11 has incorrect dimension. Defaulting to zero')
end
if ~exist('B12','var')
    B12=zeros(nx,nu);
    disp('B12 is undefined. Defaulting to zero')
elseif any(size(B12)~=[nx,nu])
    B12=zeros(nx,nu);
    disp('B12 has incorrect dimension. Defaulting to zero')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if PDE parameters are defined and are of appropriate
% dimensions
if ~exist('A0','var')
    disp('A0 is undefined. Defaulting to zero');
    A0=zeros(np);
elseif any(size(A0)~= [np,np])
    disp('A0 has incorrect dimension. Defaulting to zero');
    A0=zeros(np);
end

if ~exist('A1','var')
    disp('A1 is undefined. Defaulting to zero');
    A1=zeros(np,n1+n2);
elseif any(size(A1)~= [np,n1+n2])
    disp('A1 has incorrect dimension. Defaulting to zero');
    A1=zeros(np,n1+n2);
end

if ~exist('A2','var')
    disp('A2 is undefined. Defaulting to zero');
    A2=zeros(np,n2);
elseif any(size(A2)~=[np,n2])
    disp('A2 has incorrect dimension. Defaulting to zero');
    A2=zeros(np,n2);
end

if ~exist('E','var')
    disp('E is undefined. Defaulting to zero');
    E = zeros(np,nx);
elseif any(size(E)~=[np,nx])
    disp('E has incorrect dimension. Defaulting to zero');
    E = zeros(np,nx);
end

if ~exist('B21','var')
    B21=zeros(np,nw);
    disp('B21 is undefined. Defaulting to zero')
elseif any(size(B21)~=[np,nw])
    B21=zeros(np,nw);
    disp('B21 has incorrect dimension. Defaulting to zero')
end

if ~exist('B22','var')
    B22=zeros(np,nu);
    disp('B22 is undefined. Defaulting to zero')
elseif any(size(B22)~=[np,nu])
    B22=zeros(np,nu);
    disp('B22 has incorrect dimension. Defaulting to zero')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if BCs are defined and are of appropriate
% dimensions. 
% NOTE: B is a mandatory input
if ~exist('B','var')
    error('B is undefined. Boundary conditions are necessary for PDEs');
elseif any(size(B)~=[n1+2*n2,2*n1+4*n2])
    error('B has incorrect dimension. Not enough boundary conditions');
end

if ~exist('Bx','var')
    disp('Bx is undefined. Defaulting to zero');
    Bx = zeros(n1+2*n2,nx);
elseif any(size(Bx)~=[n1+2*n2,nx])
    disp('Bx has incorrect dimension. Defaulting to zero');
    Bx = zeros(n1+2*n2,nx);
end

if ~exist('Bw','var')
    disp('Bw is undefined. Defaulting to zero');
    Bw = zeros(n1+2*n2,nw);
elseif any(size(Bw)~=[n1+2*n2,nw])
    disp('Bw has incorrect dimension. Defaulting to zero');
    Bw = zeros(n1+2*n2,nw);
end

if ~exist('Bu','var')
    disp('Bu is undefined. Defaulting to zero');
    Bu = zeros(n1+2*n2,nu);
elseif any(size(Bu)~=[n1+2*n2,nu])
    disp('Bu has incorrect dimension. Defaulting to zero');
    Bu = zeros(n1+2*n2,nu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if regulated output parameters are defined and are of appropriate
% dimensions
if ~exist('C1','var')
    C1=zeros(nz,nx);
    disp('C1 is undefined. Defaulting to zero')
elseif any(size(C1)~=[nz,nx])
    C1=zeros(nz,nx);
    disp('C1 has incorrect dimension. Defaulting to zero')
end
if ~exist('C10','var')
    C10=zeros(nz,2*n1+4*n2);
    disp('C10 is undefined. Defaulting to zero')
elseif any(size(C10)~=[nz,2*n1+4*n2])
    C10=zeros(nz,2*n1+4*n2);
    disp('C10 has incorrect dimension. Defaulting to zero')
end
if ~exist('Ca1','var')
    Ca1=zeros(nz,np);
    disp('Ca1 is undefined. Defaulting to zero')
elseif any(size(Ca1)~=[nz,np])
    Ca1=zeros(nz,np);
    disp('Ca1 has incorrect dimension. Defaulting to zero')
end
if ~exist('Cb1','var')
    Cb1=zeros(nz,n1+n2);
    disp('Cb1 is undefined. Defaulting to zero')
elseif any(size(Cb1)~=[nz,n1+n2])
    Cb1=zeros(nz,n1+n2);
    disp('Cb1 has incorrect dimension. Defaulting to zero')
end
if ~exist('Cc1','var')
    Cc1=zeros(nz,n2);
    disp('Cc1 is undefined. Defaulting to zero')
elseif any(size(Cc1)~=[nz,n2])
    Cc1=zeros(nz,n2);
    disp('Cc1 has incorrect dimension. Defaulting to zero')
end
if ~exist('D11','var')
    D11=zeros(nz,nw);
    disp('D11 is undefined. Defaulting to zero')
elseif any(size(D11)~=[nz,nw])
    D11=zeros(nz,nw);
    disp('D11 has incorrect dimension. Defaulting to zero')
end
if ~exist('D12','var')
    D12=zeros(nz,nu);
    disp('D12 is undefined. Defaulting to zero')
elseif any(size(D12)~=[nz,nu])
    D12=zeros(nz,nu);
    disp('D12 has incorrect dimension. Defaulting to zero')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we check if sensed output parameters are defined and are of appropriate
% dimensions
if ~exist('C2','var')
    C2=zeros(ny,nx);
    disp('C2 is undefined. Defaulting to zero')
elseif any(size(C2)~=[ny,nx])
    C2=zeros(ny,nx);
    disp('C2 has incorrect dimension. Defaulting to zero')
end
if ~exist('C20','var')
    C20=zeros(ny,2*n1+4*n2);
    disp('C20 is undefined. Defaulting to zero')
elseif any(size(C20)~=[ny,2*n1+4*n2])
    C20=zeros(ny,2*n1+4*n2);
    disp('C20 has incorrect dimension. Defaulting to zero')
end
if ~exist('Ca2','var')
    Ca2=zeros(ny,np);
    disp('Ca2 is undefined. Defaulting to zero')
elseif any(size(Ca2)~=[ny,np])
    Ca2=zeros(ny,np);
    disp('Ca2 has incorrect dimension. Defaulting to zero')
end
if ~exist('Cb2','var')
    Cb2=zeros(ny,n1+n2);
    disp('Cb2 is undefined. Defaulting to zero')
elseif any(size(Cb2)~=[ny,n1+n2])
    Cb2=zeros(ny,n1+n2);
    disp('Cb2 has incorrect dimension. Defaulting to zero')
end
if ~exist('Cc2','var')
    Cc2=zeros(ny,n2);
    disp('Cc2 is undefined. Defaulting to zero')
elseif any(size(Cc2)~=[ny,n2])
    Cc2=zeros(ny,n2);
    disp('Cc2 has incorrect dimension. Defaulting to zero')
end
if ~exist('D21','var')
    D21=zeros(ny,nw);
    disp('D21 is undefined. Defaulting to zero')
elseif any(size(D21)~=[ny,nw])
    D21=zeros(ny,nw);
    disp('D21 has incorrect dimension. Defaulting to zero')
end
if ~exist('D22','var')
    D22=zeros(ny,nu);
    disp('D22 is undefined. Defaulting to zero')
elseif any(size(D22)~=[ny,nu])
    D22=zeros(ny,nu);
    disp('D22 has incorrect dimension. Defaulting to zero')
end