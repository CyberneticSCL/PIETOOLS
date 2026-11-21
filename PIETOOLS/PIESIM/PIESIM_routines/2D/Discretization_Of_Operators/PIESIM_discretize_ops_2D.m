%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_ops_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of PIE operators with Chebyshev methods  
%
% Inputs:
% 1) PIE - PIE structure of the problem
%    ('struct' with the 2D operators T, Tw, Tu, A, B1, B2, C1, C2, D11,
%    D12, D21, D22, each a 'struct' of its kernels R00, R0x, ..., R22 with
%    fields dim, var1, var2; and misc, with the optional map Tmap)
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
%    ('struct' with fields N = [N_x, N_y], nx, ny, n)
%
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% all PIE operators
%    ('struct' with fields D11cheb, D12cheb, D21cheb, D22cheb, C1cheb,
%    C2cheb, Tucheb, Twcheb, B1cheb, B2cheb, Acheb; Tcheb_inv = inv(T);
%    Atotal = inv(T)*A; Tcheb_2PDEstate, Tucheb_2PDEstate,
%    Twcheb_2PDEstate and, if PIE.misc has the field Tmap,
%    Tchebmap_2PDEstate, the maps to the primary states, not truncated)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_ops_2D
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
% Initial coding YP  - 04_16_2024
% YP 06/26/2025 - modified support for T0 operator with new pie_struct
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output fields renamed –
%   Mcheb_inv to Tcheb_inv, Mcheb_nonsquare, Tucheb_nonsquare,
%   Twcheb_nonsquare to Tcheb_2PDEstate, Tucheb_2PDEstate, Twcheb_2PDEstate;
%   the optional map is PIE.misc.Tmap, output Tchebmap_2PDEstate (was
%   PIE.misc.T0, Mcheb0_nonsquare); C1, C2 discretized with flag 1 and
%   Tu, Tw, B1, B2 with flag 2 (were 2 and 1); help and comments revised.

function Dop=PIESIM_discretize_ops_2D(PIE, psize);

% Last argument (flag): 0 maps inputs to outputs, 1 states to outputs, 2
% inputs to states, 3 states to states
 Dop.D11cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D11,psize,0);
 Dop.D12cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D12,psize,0);
 Dop.D21cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D21,psize,0);
 Dop.D22cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D22,psize,0);
 Dop.C1cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.C1,psize,1);
 Dop.C2cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.C2,psize,1);
 % Input-to-state operators; Tu, Tw also as maps to the primary states
 [Dop.Tucheb,Dop.Tucheb_2PDEstate]=PIESIM_fullPI2Mat_cheb_2D(PIE.Tu,psize,2);
 [Dop.Twcheb,Dop.Twcheb_2PDEstate]=PIESIM_fullPI2Mat_cheb_2D(PIE.Tw,psize,2);
 Dop.B1cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.B1,psize,2);
 Dop.B2cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.B2,psize,2);
 % State-to-state operators; T (and Tmap) also as maps from the
 % fundamental to the primary states
 Dop.Acheb=PIESIM_fullPI2Mat_cheb_2D(PIE.A,psize,3);
 [Tcheb, Dop.Tcheb_2PDEstate]=PIESIM_fullPI2Mat_cheb_2D(PIE.T,psize,3);
 if isfield(PIE.misc,'Tmap') 
     [Tchebmap, Dop.Tchebmap_2PDEstate]=PIESIM_fullPI2Mat_cheb_2D(PIE.misc.Tmap,psize,3);
 end
%  
  % Time propagator of the discretized PIE T*a_t = A*a: Atotal = inv(T)*A
  Dop.Tcheb_inv=inv(Tcheb);
  Dop.Atotal=Dop.Tcheb_inv*Dop.Acheb;





