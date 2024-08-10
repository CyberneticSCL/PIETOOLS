%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_ops_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of PIE operators with Chebyshev methods  
%
% Inputs:
% 1) PIE - PIE structure of the problem
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
%
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% all PIE operators
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 04_16_2024

function Dop=PIESIM_discretize_ops_2D(PIE, psize);

% Define local variables

no=psize.no;
nw=psize.nw;
nu=psize.nu;
nro=psize.nro;
noo=psize.noo;
N=psize.N;

% Define "Degree of smoothness for each PDE state". Currently we keep track of it in a discretization of PIE
% for a consistent representation of the original states. This is needed to properly
% account for the number of Chebyshev coefficients that are being solved
% for for each state: N+1 coefficients for n0 states (polynomilas of degree N), N coefficients for n1
% states (polynomials of degree N-1), N-1 coefficients for n2 states
% (polynomials of degree N-2) etc.


disp('Setting up Chebyshev matrices for the PIE system');

% Discretize 4PI operators

% Set the last entry to PIESIM_4PI2Mat_cheb to 1 if a structure is a full 4PI operator (for A and T)

% Set the last entry to PIESIM_4PI2Mat_cheb to 0 if a structure has an empty right side (for Tw, Tu, B1 and B2
% operators)

% Set the last entry to PIESIM_4PI2Mat_cheb to 2 if a structure has an
% empty bottom row (for C1 and C2 operators)

 Dop.Twcheb=PIESIM_fullPI2Mat_cheb_2D(nw,PIE.Tw,psize,0);
 Dop.Tucheb=PIESIM_fullPI2Mat_cheb_2D(nu,PIE.Tu,psize,0);
 Dop.B1cheb=PIESIM_fullPI2Mat_cheb_2D(nw,PIE.B1,psize,0);
 Dop.B2cheb=PIESIM_fullPI2Mat_cheb_2D(nu,PIE.B2,psize,0);
 Dop.C1cheb=PIESIM_fullPI2Mat_cheb_2D(nro,PIE.C1,psize,2);
 Dop.C2cheb=PIESIM_fullPI2Mat_cheb_2D(noo,PIE.C2,psize,2);
 Dop.Acheb=PIESIM_fullPI2Mat_cheb_2D(no,PIE.A,psize,1);
 [Mcheb, Dop.Mcheb_nonsquare]=PIESIM_fullPI2Mat_cheb_2D(no,PIE.T,psize,1);
 if isfield(PIE,'T0') 
     [Mcheb0, Dop.Mcheb0_nonsquare]=PIESIM_fullPI2Mat_cheb_2D(no,PIE.T0,psize,1);
 end
%  
  Dop.Mcheb_inv=inv(Mcheb);
  Dop.Atotal=Dop.Mcheb_inv*Dop.Acheb;





