%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_ops.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of PIE operators with Chebyshev methods  
%
% Inputs:
% 1) PIE - PIE structure of the problem, with PI operators, T,Tu,Tw,A,Bi currently supported
% 2) psize - size of the problem: contains the variables nu,nw,nx,N,n0,n1,n2
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for T,Tu,Tw,A,Bi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 11_05_2021

function Dop=PIESIM_discretize_ops(PIE, psize);

% Define local variables

nx=psize.nx;
nw=psize.nw;
nu=psize.nu;
N=psize.N;

% Define "Degree of smoothness for each PDE state". Currently we keep track of it in a discretization of PIE
% for a consistent representation of the original states. This is needed to properly
% account for the number of Chebyshev coefficients that are being solved
% for for each state: N+1 coefficients for n0 states (polynomilas of degree N), N coefficients for n1
% states (polynomials of degree N-1), N-1 coefficients for n2 states
% (polynomials of degree N-2) etc.

% Define degree of smoothness p

psize_aux1=[1 psize.n];
psize_aux0=[0 psize.n];
nsum=cumsum(psize.n);
nsump1=cumsum(psize_aux1);
nsump0=cumsum(psize_aux0);
for i=1:length(psize.n);
p(nsump1(i):nsum(i))=i-1;
end


disp('Setting up Chebyshev matrices for the PIE system');

% Discretize 4PI operators

% Set the last entry to PIESIM_4PI2Mat_cheb to 1 if a structure is a full 4PI operator (for A and T)

% Set the last entry to PIESIM_4PI2Mat_cheb to 0 if a structure has an empty right side (for Tw, Tu, B1 and B2
% operators)

 Dop.Twcheb=PIESIM_4PI2Mat_cheb(N,nw,PIE.Tw,p,0);
 Dop.Tucheb=PIESIM_4PI2Mat_cheb(N,nu,PIE.Tu,p,0);
 Dop.B1cheb=PIESIM_4PI2Mat_cheb(N,nw,PIE.B1,p,0);
 Dop.B2cheb=PIESIM_4PI2Mat_cheb(N,nu,PIE.B2,p,0);
 Dop.Acheb=PIESIM_4PI2Mat_cheb(N,nx,PIE.A,p,1);
[Mcheb, Dop.Mcheb_nonsquare]=PIESIM_4PI2Mat_cheb(N,nx,PIE.T,p,1);
%  
  Dop.Mcheb_inv=inv(Mcheb);
  Dop.Atotal=Dop.Mcheb_inv*Dop.Acheb;




