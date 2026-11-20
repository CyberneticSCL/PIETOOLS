%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_ops.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of PIE operators with Chebyshev methods  
%
% Inputs:
% 1) PIE - PIE structure of the problem
%    ('struct' with the operators T, Tw, Tu, A, B1, B2, C1, C2, D11, D12,
%    D21, D22, each a 'struct' of its kernels P, Q1, Q2, R.R0, R.R1, R.R2
%    with the field dim; and misc, with the optional maps Tmap, Tumap,
%    Twmap)
% 2) psize - size of the problem: contains the variables N, n
%    ('struct'; N the Chebyshev order, n(i) the number of PDE states
%    differentiable i-1 times)
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% T,Tu,Tw,A,Bi,Ci, Dij
%    ('struct' with fields D11cheb, D12cheb, D21cheb, D22cheb, C1cheb,
%    C2cheb, Tucheb, Twcheb, B1cheb, B2cheb, Acheb; Tcheb_inv = inv(T);
%    Atotal = inv(T)*A; Tcheb_2PDEstate, Tucheb_2PDEstate,
%    Twcheb_2PDEstate, the maps to the primary states, not truncated)
% Optional operators (when PIE.misc has the field Tmap, Tumap or Twmap) are
% Dop.Tchebmap_2PDEstate, Dop.Tuchebmap_2PDEstate and Dop.Twchebmap_2PDEstate
% In this case, T, Tu, Tw signify LHS dynamics, and Tmap, Tumap, Twmap are the
% mapping operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_ops
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
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
% MMP, 11/19/2026: Release comment pass (comment-only).

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


 % Input-to-output and state-to-output operators
 Dop.D11cheb=PIESIM_4PI2Mat_cheb(N,PIE.D11,p,0);
 Dop.D12cheb=PIESIM_4PI2Mat_cheb(N,PIE.D12,p,0);
 Dop.D21cheb=PIESIM_4PI2Mat_cheb(N,PIE.D21,p,0);
 Dop.D22cheb=PIESIM_4PI2Mat_cheb(N,PIE.D22,p,0);
 Dop.C1cheb=PIESIM_4PI2Mat_cheb(N,PIE.C1,p,1);
 Dop.C2cheb=PIESIM_4PI2Mat_cheb(N,PIE.C2,p,1);
 % Input-to-state operators; Tu, Tw (and Tumap, Twmap) also as maps to the
 % primary states
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
 % State-to-state operators; T (and Tmap) also as maps from the
 % fundamental to the primary states
 Dop.Acheb=PIESIM_4PI2Mat_cheb(N,PIE.A,p,3);
 [Tcheb, Dop.Tcheb_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.T,p,3);
 if isfield(PIE.misc,'Tmap') 
     [Tchebmap, Dop.Tchebmap_2PDEstate]=PIESIM_4PI2Mat_cheb(N,PIE.misc.Tmap,p,3);
 end
%  
  % Time propagator of the discretized PIE T*a_t = A*a: Atotal = inv(T)*A
  Dop.Tcheb_inv=inv(Tcheb);
  Dop.Atotal=Dop.Tcheb_inv*Dop.Acheb;




