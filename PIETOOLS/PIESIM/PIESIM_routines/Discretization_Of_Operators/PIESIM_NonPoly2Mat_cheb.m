%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_NonPoly2Mat_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a non-polymonial operator
% Discretizes FULL OPERATOR but does not provide a PIE-2-PDE state map 
%
% Inputs:
% 1) N   - polynomial order of Chebyshev discretization polynomial
% 2) Rop -  non-polynomial matrix operator 
% 3) p - scalar - a "degree of smoothness" vector
% 4) grid - field containing the following sub-fields:
% grid.phys - physical grid for states differentiable up to order zero (corresponding to a primary = PDE state discretization)
% grid - cell array containing grids of different degrees of differentiability
%
% Outputs:
% 1) A - matrix that represents a Chebyshev
% discretizaiton of a non-polynomial matrix operator. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_NonPoly2Mat_cheb
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
% Initial coding YP  - 11_27_2021
% DJ, 12/16/2024: Remove redundant variables.
% MMP, 11/19/2026: Release comment pass (comment-only).

function A=PIESIM_NonPoly2Mat_cheb(N, Rop, p, grid)

ns=size(p,2);
no=size(Rop,2);

for m=ns:-1:1
    rsize=N-p(m)+1;

    Alocal(1:rsize,1:no)=0;
    for i=1:no 
        var_force=double(subs(Rop(m,i),grid{p(m)+1}));
        Alocal(:,i)=fcht(var_force);
    end % i loop

    A_cell{m,1}=Alocal;

end % m loop

A=cell2mat(A_cell);
