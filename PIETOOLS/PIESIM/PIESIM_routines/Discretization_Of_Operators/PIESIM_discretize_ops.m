%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_ops.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of PIE operators with Chebyshev methods  
%
% Inputs:
% 1) PIE - PIE structure of the problem
% 2) psize - size of the problem: contains the variables nu,nw,no,N,n0,n1,n2
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% T,Tu,Tw,A,Bi,Ci, Dij
% Optional operators (when PIE.misc exists) are Dop.Tmap, Dop.Tumap and Dop. Twmap (all with cheb extensions) 
% In this case, T, Tu, Tw signify LHS dynamics, and Tmap, Tumap, Twmap are the
% mapping operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 11_05_2021
% YP - 6/16/2022 - added discretization of C1, C2 operators to allow for
% computation of observed and regulated outputs. Separated T (PIE operator that enters the LHS of the dynamics) 
% from Tmap (fundamental-state to PDE-state mapping operator). Needed if the operator
% in the LHS of the time propagator is different from the mapping operator
% (e.g. in weighted formulation for cylindrical coordinates when PDE is multipled by a power of r).
% YP - added functionality to support infinite-dimensional disturbances
% through parser - 6_1_2025
% YP 06/26/2025 - modified support for T0 operator with new pie_struct
% YP 1/6/2026 - modified Tu, Tw operators to allow reconstruction of disturbance-to-state
% maps, added support for Tumap and Twmap operators (needed for cylindrical coordinates) 

function Dop=PIESIM_discretize_ops(PIE, psize);

% Define local variables

N=psize.N;

% Define "Degree of smoothness for each PDE state". Currently we keep track of it in a discretization of PIE
% for a consistent representation of the original states. This is needed to properly
% account for the number of Chebyshev coefficients that are being solved
% for for each state: N+1 coefficients for n0 states (polynomilas of degree N), N coefficients for n1
% states (polynomials of degree N-1), N-1 coefficients for n2 states
% (polynomials of degree N-2) etc.

% Define degree of smoothness p

p = repelem(0:length(psize.n)-1, psize.n);

%disp('Setting up Chebyshev matrices for the PIE system');

% Discretize 4PI operators

 % The rightmost entry to PIESIM_4PI2Mat_cheb is a flag.
% flag = 0 if a structure maps disturbances or control inputs to regulated or observed outputs (for D11, D12,
% D21, D22 operators)
% flag = 1 if a structure maps solution states (ODE+PDE) to regulated or observed outputs (for C1, C2 operators)
% flag = 2 if a structure maps disturbances or control inputs to solution states (ODE+PDE) (for Tw, Tu, B1 and B2
% operators)
% flag = 3 if a structure maps solution states (ODE+PDE) to solution states (ODE+PDE) (for A and T)
%


 Dop.D11cheb=PIESIM_4PI2Mat_cheb(N,PIE.D11,p,0);
 Dop.D12cheb=PIESIM_4PI2Mat_cheb(N,PIE.D12,p,0);
 Dop.D21cheb=PIESIM_4PI2Mat_cheb(N,PIE.D21,p,0);
 Dop.D22cheb=PIESIM_4PI2Mat_cheb(N,PIE.D22,p,0);
 Dop.C1cheb=PIESIM_4PI2Mat_cheb(N,PIE.C1,p,1);
 Dop.C2cheb=PIESIM_4PI2Mat_cheb(N,PIE.C2,p,1);
 [Dop.Tucheb,Dop.Tucheb_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.Tu,p,2);
 [Dop.Twcheb,Dop.Twcheb_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.Tw,p,2);
 if isfield(PIE.misc,'Tumap') 
     [Tuchebmap,Dop.Tuchebmap_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.misc.Tumap,p,2);
 end
 if isfield(PIE.misc,'Twmap') 
     [Twchebmap,Dop.Twchebmap_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.misc.Twmap,p,2);
 end
 Dop.B1cheb=PIESIM_4PI2Mat_cheb(N,PIE.B1,p,2);
 Dop.B2cheb=PIESIM_4PI2Mat_cheb(N,PIE.B2,p,2);
 Dop.Acheb=PIESIM_4PI2Mat_cheb(N,PIE.A,p,3);
 [Tcheb, Dop.Tcheb_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.T,p,3);
 if isfield(PIE.misc,'Tmap') 
     [Tchebmap, Dop.Tchebmap_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.misc.Tmap,p,3);
 end
%  
  Dop.Tcheb_inv=inv(Tcheb);
  Dop.Atotal=Dop.Tcheb_inv*Dop.Acheb;




