%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_4PI2Mat_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs a discrete version of the multi-dimensional 4PI operator 
% Discretizes FULL ASSEMLED OPERATOR
%   [x0; x1] -> [P x0 + int_{-1}^{1} Q1(t) x1(t) dt;
%                Q2(s) x0 + R0(s) x1(s) + int_{-1}^{s} R1(s,t) x1(t) dt
%                         + int_{s}^{1} R2(s,t) x1(t) dt],   s in [-1,1],
% as A = [P Q1; Q2 R] acting on the values of x0 and the Chebyshev
% coefficients of x1.
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop - 4PI operator 
%       ('struct' with fields P, Q1, Q2 and R (fields R0, R1, R2))
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for 3PI: p(j) is the degree of
% differentiability of state j, whose Chebyshev expansion has degree N-p(j)
%
% flag = 0 if a structure maps disturbances or control inputs to regulated or observed outputs (for D11, D12,
% D21, D22 operators)
% flag = 1 if a structure maps solution states (ODE+PDE) to regulated or observed outputs (for C1, C2 operators)
% flag = 2 if a structure maps disturbances or control inputs to solution states (ODE+PDE) (for Tw, Tu, B1 and B2
% operators)
% flag = 3 if a structure maps solution states (ODE+PDE) to solution states (ODE+PDE) (for A and T)
%
%
%
% Outputs:
% A - discretization matrix of the time-dependent PIE propagator
% A_2PDEstate -discretization matrix for transform between
% fundamental and primary solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_4PI2Mat_cheb
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
% Initial coding YP  - 5_31_2021
% YP - added functionality to support infinite-dimensional disturbances and
% regulated/observed outputs - 6_2_2025
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
%   renamed A_2PDEstate; help and comments revised.

function [A, A_2PDEstate]=PIESIM_4PI2Mat_cheb(N, Rop, p, flag)

% Iniatilize the matrix blocks of the full PI operator as empty
 Pblock=[];
 Q1block=[];
 Q2block=[];
 Rblock=[];

 Q2block_2PDEstate=[];
 Rblock_2PDEstate=[];

% Pblock: size nscalar1 x nscalar2
Pblock=double(Rop.P);

% Q1block: size nscalar1 x nspat2


if (size(Rop.Q1)~=0)

if(mod(flag,2)==0)
pcol=zeros(size(Rop.Q1,2),1);
else
pcol=p;
end

Q1block=PIESIM_PI2Mat_opint_cheb(N, Rop.Q1, pcol);
end


% Q2block: size nspat1 x nscalar2

% Q2block_2PDEstate: nspat1_full x nscalar2


if (size(Rop.Q2)~=0)

    if (flag<=1)
        prow=zeros(size(Rop.Q2,1),1);
    else
        prow=p;
    end


if (nargout==1)
 Q2block= PIESIM_Poly2Mat_cheb(N, Rop.Q2, prow);
else
[Q2block, Q2block_2PDEstate]= PIESIM_Poly2Mat_cheb(N, Rop.Q2, prow);
end

end % size(Rop.Q2)~=0


% Rblock: size  nspat1 x nspat2


if (size(Rop.R.R0)~=0)

% Smoothness of the rows (states if flag >= 2) and of the columns (states
% if flag is odd); zero for outputs and inputs
switch flag
    case 0
        prow=zeros(size(Rop.R.R0,1),1);
        pcol=zeros(size(Rop.R.R0,2),1);
    case 1
        prow=zeros(size(Rop.R.R0,1),1);
        pcol=p;
    case 2
        prow=p;
        pcol=zeros(size(Rop.R.R0,2),1);
    case 3
        prow=p;
        pcol=p;

end


% 3PI block, and its map to the primary states if requested
if (nargout==1)
Rblock=PIESIM_3PI2Mat_cheb(N, Rop.R, prow, pcol);
else
[Rblock, Rblock_2PDEstate]=PIESIM_3PI2Mat_cheb(N, Rop.R, prow, pcol);
end

end % size(Rop.R.R0,2)~=0


A=[Pblock Q1block; Q2block Rblock];
if (nargout>1)
    A_2PDEstate=[Pblock Q1block; Q2block_2PDEstate Rblock_2PDEstate];
end





