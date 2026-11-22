%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Poly2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a polymonial operator in
% 2D
% Discretizes a FULL OPERATOR
% A maps v in R^no to the Chebyshev coefficients of R(x,y) v, x,y in
% [-1,1]: column k holds, for each row m of R, the coefficients a_ij of
% R(m,k)(x,y) = sum_{i<=N_x, j<=N_y} a_ij T_i(x) T_j(y), x index i fastest.
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
%       (1x2 vector [N_x, N_y], the orders in x and y)
% Rop -  polynomial matrix operator 
%       (ns x no 'polynomial' or 'double' matrix R(x,y))
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% p - 2 x ns array - a "degree of smoothness" vector [p_x; p_y] per row
%     state (row of R)
%
% Outputs:
% A - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% total matrix operator for PIE solution
%     (no columns; for row state m, the rows with i <= N_x - p(1,m) and
%     j <= N_y - p(2,m))
% A_2PDEstate - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% 2PDEstate total matrix operator for reconstruction of the primary (PDE)
% solution
%     (ns*prod(N+1) x no, not truncated)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_Poly2Mat_cheb_2D
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
% DJ, 07/17/2024: bugfix in case "Rop" is not polynomial;
% DJ, 12/16/2024: Correct limits of for loop j=1:size(acheb);
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.
% YP, MMP, 11/22/2026: Changes since PIETOOLS 2025: output A_nonsquare
%   renamed A_2PDEstate (ns*prod(N+1) rows, not truncated); N is [N_x, N_y]
%   and p is a 2 x ns array [p_x; p_y]; R is evaluated on deg+2
%   Chebyshev-Gauss-Lobatto nodes per direction (deg the maximal degree of
%   R), its coefficients padded or truncated to (N_x+1) x (N_y+1); A is
%   A_2PDEstate with the rows of each state truncated by p.

function [A, A_2PDEstate]=PIESIM_Poly2Mat_cheb_2D(N, Rop, var1, p)

% Extract spatial variables (x,y) defining the operator R                   % DJ, 12/16/2024
s1 = var1(1);   s2 = var1(2);

% Evaluate everything on maximum grid and then compress

% ns row states, no columns of R
ns=size(p,2);
no=size(Rop,2);

% deg+2 Chebyshev-Gauss-Lobatto nodes per direction, exact for the
% maximal degree deg of R
if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));

% Block m (row state m): column k holds the coefficients of R(m,k)
for m=1:ns
Alocal(1:prod(N+1),1:no)=0;

for k=1:no
    % Values of R(m,k) on the 2D node grid, rows x, columns y
    if (isempty(Rop))
        Rop_local=0;
    else
        Rop_local=Rop(m,k);
    end
    if isa(Rop_local,'polynomial')
    Reval=subs(subs(Rop_local,s1,chebgrid)',s2,chebgrid);
    else
    Reval=Rop_local*ones(deg+2,deg+2);
    end

    % Initiailize square array of size (N(1)+1)x(N(2)+1)
    acheb2D_all=zeros(N(1)+1,N(2)+1);
    acheb2D=fcgltran2d(double(Reval),1);
    % Coefficients of orders up to N_x, N_y; higher orders, present when    % YP, MMP, 11/22/2026
    % deg >= N, are dropped, so the array keeps its size                    % YP, MMP, 11/22/2026
    nx = min(deg+2,N(1)+1);   ny = min(deg+2,N(2)+1);                       % YP, MMP, 11/22/2026
    acheb2D_all(1:nx,1:ny) = acheb2D(1:nx,1:ny);                            % YP, MMP, 11/22/2026

    % Coefficients a_ij, padded to (N_x+1)x(N_y+1), as column k (x index
    % fastest)
    acheb=reshape(acheb2D_all,[],1);

for j=1:length(acheb)                                                     
Alocal(j,k)=acheb(j);
end

end % k loop

% Truncate time propagator matrix to the correct dimensions
A_propagator=PIESIM_Mat_Truncate(Alocal,N,p(:,m),'row');
A_cell{m,1}=A_propagator;
A_cell_2PDEstate{m,1}=Alocal;

end % m loop


A_2PDEstate=cell2mat(A_cell_2PDEstate);
A=cell2mat(A_cell);
