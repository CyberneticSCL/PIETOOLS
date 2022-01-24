%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PDE_batch2terms.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDE_terms=convert_PIETOOLS_PDE_batch2terms(PDE_batch)
% converts the PDE object format from batch to terms. Refer to files
% initialize_PIETOOLS_PDE_batch and initialize_PIETOOLS_PDE_terms for more
% information on those input formats.
PDE=initialize_PIETOOLS_PDE_batch(PDE_batch);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-Defined Inputs for specifying the dynamics:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: All dimensions (nx,nw,nu,nz,ny,n0,n1,n2) are assumed to be 0 by default.
% Specifically, (all optional)
% PDE.nx     -  number of ODE states
% PDE.nw     -  number of disturbances
% PDE.nu     -  number of controlled inputs
% PDE.nz     -  number of regulated outputs
% PDE.ny     -  number of observed outputs
% PDE.n0     -  number of undifferentiated PDE states
% PDE.n1     -  number of once spatially differentitaed PDE states
% PDE.n2     -  number of twice spatially differentiated PDE states
%
% NOTE: Any matrix with a 0-dimension should be ommitted
%
% PDE Terms (required)
% PDE.dom    -  interval of the domain of the spatial variable - s \in [a,b] (Required)

% PDE.A  - nx x nx matrix  
% PDE.E0 - nx x (2*n1+4*n2) matrix
% PDE.Ea - nx x (n0+n1+n2) matrix valued polynomial in s
% PDE.Eb - nx x (n1+n2) matrix valued polynomial in s
% PDE.Ec - nx x (n2) matrix valued polynomial in s
% PDE.B11 - nx x nw matrix 
% PDE. B12 - nx x nu matrix
% 
% PDE.E - (n0+n1+n2) x nx matrix valued polynomial in s
% PDE.A0 - (n0+n1+n2) x (n0+n1+n2) matrix valued polynomial in s
% PDE.A1 - (n0+n1+n2) x (n1+n2) matrix valued polynomial in s
% PDE.A2 - (n0+n1+n2) x (n2) matrix valued polynomial in s
% PDE.B21 - (n0+n1+n2) x nw matrix valued polynomial in s
% PDE.B22 - (n0+n1+n2) x nu matrix valued polynomial in s
% 
% PDE.B - (n1+2*n2) x (2*n1+4*n2) matrix
% PDE.Bx - (n1+2*n2) x nx matrix
% PDE.Bw - (n1+2*n2) x nw matrix
% PDE.Bu - (n1+2*n2) x nu matrix
% 
% PDE.C1 - nz x nx matrix
% PDE.C2 - ny x nx matrix
% PDE.C10 - nz x (2*n1+4*n2) matrix
% PDE.C20 - ny x (2*n1+4*n2) matrix
% PDE.Ca1 - nz x (n0+n1+n2) matrix
% PDE.Ca2 - ny x (n0+n1+n2) matrix
% PDE.Cb1 - nz x (n1+n2) matrix
% PDE.Cb2 - ny x (n1+n2) matrix
% PDE.Cc1 - nz x n2 matrix
% PDE.Cc2 - ny x n2 matrix
% PDE.D11 - nz x nw matrix
% PDE.D12 - nz x nu matrix
% PDE.D21 - ny x nw matrix
% PDE.D22 - ny x nu matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS  - 5_30_2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for domain (mandatory variables)
if ~isfield(PDE,'dom')
    PDE.dom = [0,1];
    disp('Warning: PDE domain not defined. Defaulting to [0,1]');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize proper batch PDE structure if not initialized
PDE = initialize_PIETOOLS_PDE_batch(PDE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we initialize a terms structure PDE with correct dimensions and
% domain.
PDE_out.n.nx = PDE.nx;
PDE_out.n.n_pde = [PDE.n0 PDE.n1 PDE.n2];
PDE_out.n.nz = PDE.nz;
PDE_out.n.ny = PDE.ny;
PDE_out.n.nw = PDE.nw;
PDE_out.n.nu = PDE.nu;
PDE_out.n.nv = PDE.nx+PDE.nw+PDE.nu; 
PDE_out.n.nr = PDE.nx+PDE.nz+PDE.ny + 2*PDE.n1 + 4*PDE.n2; %?

PDE_out.dom = PDE.dom;

PDE_out = initialize_PIETOOLS_PDE_terms(PDE_out);

N = length(PDE_out.n.n_pde)-1;
% find total PDE states
np = sum(PDE_out.n.n_pde(1:N+1));
% find all possible derivatives length
% % np_all_derivatives = sum((1:N+1).*PDE.n.n_pde);
np_all_derivatives = sum(1:N+1);
% find total possible Boundary values
nBVs=2*sum((0:N).*PDE_out.n.n_pde); nBC = nBVs/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define number of variable explicitly for further reference
nx = PDE_out.n.nx;  n0 = PDE.n0; n2 = PDE.n2; n1 = PDE.n1;
nw = PDE_out.n.nw;  nu = PDE_out.n.nu;  nr = PDE_out.n.nr;
nz = PDE_out.n.nz;  ny = PDE_out.n.ny;  nv = PDE_out.n.nv;

% Assuming dimensions are now known, initialize the remaining parameters to zero based on
% the dimensions
pvar s theta;

PDE.vars = [s;theta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First place ODE related parameters in correct places

PDE_out.ODE.A = PDE.A;
PDE_out.ODE.Bxw = PDE.B11;
PDE_out.ODE.Bxu = PDE.B12;
PDE_out.ODE.Bxr = [PDE.E0 eye(nx) zeros(nx,nz) zeros(nx,ny)];
PDE_out.ODE.Cz = PDE.C1;
PDE_out.ODE.Dzw = PDE.D11;
PDE_out.ODE.Dzu = PDE.D12;
PDE_out.ODE.Dzr = [PDE.C10 zeros(nz,nx) eye(nz) zeros(nz,ny)];
PDE_out.ODE.Cy = PDE.C2;
PDE_out.ODE.Dyw = PDE.D21;
PDE_out.ODE.Dyu = PDE.D22;
PDE_out.ODE.Dyr = [PDE.C20 zeros(ny,nx) zeros(ny,nz) eye(ny)];
PDE_out.ODE.Cv = [eye(nx); zeros(nw,nx); zeros(nu,nx)];
PDE_out.ODE.Dvw = [zeros(nx,nw); eye(nw); zeros(nu,nw)];
PDE_out.ODE.Dvu = [zeros(nx,nu); zeros(nw,nu); eye(nu)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Place BC related parameters in correct places
% Ebp is not present in batch input, so no change needed
PDE_out.BC.Ebv = [PDE.Bx PDE.Bw PDE.Bu];

% split B*[x1(a) x1(b) x2(a) x2(b) x2s(a) x2s(b)] in Ebb structure
% j delta 0, k D = 0, l Rstate [x1 x2]
j = 0; k = 0; l = 1;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l; 
PDE_out.BC.Ebb{loc}.coeff = PDE.B(:,1:n1);
PDE_out.BC.Ebb{loc+1}.coeff = PDE.B(:,2*n1+1:2*n1+n2);
% j delta 0, k D = 1, l Rstate [x2]
j = 0; k = 1; l = 2;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l; 
PDE_out.BC.Ebb{loc}.coeff = PDE.B(:,2*n1+2*n2+1:2*n1+3*n2);
% j delta 1, k D = 0, l Rstate [x1 x2]
j = 1; k = 0; l = 1;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l; 
PDE_out.BC.Ebb{loc}.coeff = PDE.B(:,n1+1:2*n1);
PDE_out.BC.Ebb{loc+1}.coeff = PDE.B(:,2*n1+n2+1:2*n1+2*n2);
% j delta 1, k D = 1, l Rstate [x2]
j = 1; k = 1; l = 2;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l; 
PDE_out.BC.Ebb{loc}.coeff = PDE.B(:,2*n1+3*n2+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Place PDE related parameters in correct places
% Bpb is not present in batch input, so no change needed
PDE_out.PDE.Bpv = [PDE.E PDE.B21 PDE.B22];

%%% there are only 18 terms in A term structure
% j Lstate x0, l Rstate [x0 x1 x2], i I =0, k D = 0
i = 0; j = 0; k = 0; l = 0;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A0(1:n0,1:n0);
PDE_out.PDE.A{loc+1}.coeff = PDE.A0(1:n0,n0+1:n0+n1);
PDE_out.PDE.A{loc+2}.coeff = PDE.A0(1:n0,n0+n1+1:end);
% j Lstate x0, l Rstate [x1 x2], i I =0, k D = 1
i = 0; j = 0; k = 1; l = 1;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A1(1:n0,1:n1);
PDE_out.PDE.A{loc+1}.coeff = PDE.A1(1:n0,n1+1:end);
% j Lstate x0, l Rstate [x2], i I =0, k D = 2
i = 0; j = 0; k = 2; l = 2;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A2(1:n0,1:end);

%%% repeat for x1 and x2
% j Lstate x1, l Rstate [x0 x1 x2], i I =0, k D = 0
i = 0; j = 1; k = 0; l = 0;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A0(n0+1:n0+n1,1:n0);
PDE_out.PDE.A{loc+1}.coeff = PDE.A0(n0+1:n0+n1,n0+1:n0+n1);
PDE_out.PDE.A{loc+2}.coeff = PDE.A0(n0+1:n0+n1,n0+n1+1:end);
% j Lstate x1, l Rstate [x1 x2], i I =0, k D = 1
i = 0; j = 1; k = 1; l = 1;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A1(n0+1:n0+n1,1:n1);
PDE_out.PDE.A{loc+1}.coeff = PDE.A1(n0+1:n0+n1,n1+1:end);
% j Lstate x1, l Rstate [x2], i I =0, k D = 2
i = 0; j = 1; k = 2; l = 2;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A2(n0+1:n0+n1,1:end);

% j Lstate x2, l Rstate [x0 x1 x2], i I =0, k D = 0
i = 0; j = 2; k = 0; l = 0;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A0(n0+n1+1:end,1:n0);
PDE_out.PDE.A{loc+1}.coeff = PDE.A0(n0+n1+1:end,n0+1:n0+n1);
PDE_out.PDE.A{loc+2}.coeff = PDE.A0(n0+n1+1:end,n0+n1+1:end);
% j Lstate x2, l Rstate [x1 x2], i I =0, k D = 1
i = 0; j = 2; k = 1; l = 1;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A1(n0+n1+1:end,1:n1);
PDE_out.PDE.A{loc+1}.coeff = PDE.A1(n0+n1+1:end,n1+1:end);
% j Lstate x2, l Rstate [x2], i I =0, k D = 2
i = 0; j = 2; k = 2; l = 2;
loc = (i)*(N+1)*np_all_derivatives + (j)*np_all_derivatives + sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.A{loc}.coeff = PDE.A2(n0+n1+1:end,1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set r(t) = [x_b(t); int([Ea Eb Ec]*x_all_derivatives); ; int([Ca1 Cb1
% Cc1]*x_all_derivatives) ; int([Ca2 Cb2 Cc2]*x_all_derivatives)];
% PDE_out.PDE.Crp = [zeros(nBVs,np_all_ders); [Ea Eb Ec; Ca1 Cb1 Cc1; Ca2 Cb2 Cc2]];

% l Rstate [x0 x1 x2], k D = 0
k = 0; l = 0;
loc = sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.Crp{loc}.coeff = [zeros(nBVs,n0); PDE.Ea(:,1:n0);PDE.Ca1(:,1:n0);PDE.Ca2(:,1:n0)];
PDE_out.PDE.Crp{loc+1}.coeff = [zeros(nBVs,n1); PDE.Ea(:,n0+1:n0+n1);PDE.Ca1(:,n0+1:n0+n1);PDE.Ca2(:,n0+1:n0+n1)];
PDE_out.PDE.Crp{loc+2}.coeff = [zeros(nBVs,n2); PDE.Ea(:,n0+n1+1:end);PDE.Ca1(:,n0+n1+1:end);PDE.Ca2(:,n0+n1+1:end)];
% l Rstate [x1 x2], k D = 1
k = 1; l = 1;
loc = sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.Crp{loc}.coeff = [zeros(nBVs,n1);PDE.Eb(:,1:n1);PDE.Cb1(:,1:n1);PDE.Cb2(:,1:n1)];
PDE_out.PDE.Crp{loc+1}.coeff = [zeros(nBVs,n2);PDE.Eb(:,n1+1:end);PDE.Cb1(:,n1+1:end);PDE.Cb2(:,n1+1:end)];
% l Rstate [x2], k D = 2
k = 2; l = 2;
loc = sum(N:-1:N-k+1) + l + 1; 
PDE_out.PDE.Crp{loc}.coeff = [zeros(nBVs,n2);PDE.Ec(:,1:end);PDE.Cc1(:,1:end);PDE.Cc2(:,1:end)];


% Drb x_b = x_b, so set coeffs in Drb = I
% j delta 0, l Rstate [x1 x2], k D = 0
j = 0; k = 0; l = 1;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l;
PDE_out.PDE.Drb{loc}.coeff = [eye(n1); zeros(nr-n1,n1)];
PDE_out.PDE.Drb{loc+1}.coeff = [zeros(n1,n2); eye(n2); zeros(nr-n1-n2,n2)];
% j delta 0, l Rstate [x2], k D = 1
j = 0; k = 1; l = 2;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l;
PDE_out.PDE.Drb{loc}.coeff = [zeros(n1+n2,n2); eye(n2); zeros(nr-n1-2*n2,n2)];

% j delta 1, l Rstate [x1 x2], k D = 0
j = 1; k = 0; l = 1;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l;
PDE_out.PDE.Drb{loc}.coeff = [zeros(nBC,n1); eye(n1);zeros(nr-n1-nBC,n1)];
PDE_out.PDE.Drb{loc+1}.coeff = [zeros(nBC+n1,n2); eye(n2); zeros(nr-n1-n2-nBC,n2)];
% j delta 1, l Rstate [x2], k D = 1
j = 1; k = 1; l = 2;
loc = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + l;
PDE_out.PDE.Drb{loc}.coeff = [zeros(nBC+n1+n2,n2); eye(n2); zeros(nr-n1-2*n2-nBC,n2)];

PDE_out.PDE.Drv = zeros(nr,nv);
PDE_terms=PDE_out;
end 

