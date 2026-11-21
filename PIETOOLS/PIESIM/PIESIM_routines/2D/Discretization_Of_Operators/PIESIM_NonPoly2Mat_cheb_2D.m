%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_NonPoly2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a non-polymonial
% operator in 2D
%
% Inputs:
% 1) N   - polynomial order of Chebyshev discretization polynomial
% 2) Rop -  non-polynomial matrix operator 
% 3) p - scalar (non-polynomial) - a "degree of smoothness" vector
% 4) gridall - cell array of size dmax containing physical grid for all states
% depending on their degree of differentiability; dmax corresponds to the
% maximum degree of differentiability among all the states
%  gridall.x - grids in x direction
%  gridall.y - grids in y direction
%
%
% Outputs:
% A - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% square total matrix operator for PIE solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_NonPoly2Mat_cheb_2D
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
% Initial coding YP  - 4_16_2024
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: N is [N_x, N_y], each
%   state block has prod(N-p+1) rows; the function line carries the file
%   name (was PIESIM_NonPOLy2Mat_cheb_2D); help and comments revised.
function A=PIESIM_NonPoly2Mat_cheb_2D(N, Rop, p, gridall)
syms sx sy
ns=size(p,2);
no=size(Rop,2);

for m=ns:-1:1

    Alocal(1:prod(N-p(m)+1),1:no)=0;

    for i=1:no 
    var_force=double(subs(subs(Rop(m,i),sx,gridall.x{p(m)+1}),sy,gridall.y{p(m)+1}'));
    Alocal(:,i)=reshape(fcgltran2d(var_force,1),[],1);
    end % i loop

    A_cell{m,1}=Alocal;

end % m loop

A=cell2mat(A_cell);
