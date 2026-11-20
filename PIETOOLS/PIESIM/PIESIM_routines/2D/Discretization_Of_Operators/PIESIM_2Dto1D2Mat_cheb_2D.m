%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_2Dto1D2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_2PDEstate - discrete matrix representations of 2D to 1D
% operators (Rx2, Ry2)
%
%   (R v)(x) = int_{-1}^{1} [ R0(x,y) v(x,y)
%                + int_{-1}^{x} R1(x,theta,y) v(theta,y) dtheta
%                + int_{x}^{1} R2(x,theta,y) v(theta,y) dtheta ] dy
%
% for Rx2 (dir = 'x'), and likewise with x and y exchanged for Ry2. Block
% (i,j) of A_2PDEstate keeps the columns of the coefficients resolved for
% the 2D column state j; block (i,j) of A keeps, of its rows, the first
% N_s-prow(i)+1 (s = dir).
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% Rop -  a 2D to 1D component of the full 2D PI operator, corresponding to
% either Rx2 or Ry2 ('cell' {R0, R1, R2}: multiplier and the two integral
% kernels)
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% prow - vector with the degrees of differentiability for the 1D row
%        (output) states (differentiable in x for Rx2 and in y for Ry2)
% pcol - 2 x ns_col matrix with the degrees of differentiability in x and
%        y of the 2D column (input) states
% dir - direction ('x' for Rx2 and 'y' for Ry2)
%
% Outputs:
% A - Chebyshev
% discretizaiton of a Rx2/Ry2 block that represents a block of a
% total matrix operator for time-advancement of the
% spatially-discretized PIE solution 
% A_2PDEstate - Chebyshev
% discretizaiton of a Rx2/Ry2 block that represents a block of a
% total matrix operator for reconstruction of the primary (PDE)
% solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_2Dto1D2Mat_cheb_2D
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

function [A, A_2PDEstate]=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop, var1, prow, pcol, dir)

% Numbers of row and column states
ns_row=size(prow,2);
ns_col=size(pcol,2);

% Multiplier R0 and integral kernels R1, R2
  R0=Rop{1};
  R1=Rop{2};
  R2=Rop{3};

% Order of the 1D row states
if dir=='x'
    Nrow=N(1);
else
    Nrow=N(2);
end


%         % Go over the rows of the PI operator row by row and column by column

  for i=1:ns_row
     for j=1:ns_col       

% Treatment of mutliplicative operator (opmult stands for multiplicative)

   if (~isempty(R0)) 
       A_opmult=PIESIM_3PI_op_line_int_discretize_2D(N, R0(i,j), var1,dir, 0);
   end

% Treatment of integrative operators 

   if (~isempty(R1)) 
       A_opint_block_1=PIESIM_3PI_op_line_int_discretize_2D(N, R1(i,j), var1, dir, [1 0]);
   end

   if (~isempty(R2)) 
       A_opint_block_2=PIESIM_3PI_op_line_int_discretize_2D(N, R2(i,j), var1, dir, [0 1]);
   end

% Sum of the three components
   A_sum=A_opmult+A_opint_block_1+A_opint_block_2;

% Truncate the total operator to the correct dimensions

% First, PIE-to-PDE state matrix
A_2PDEstate_local=PIESIM_Mat_Truncate(A_sum,N,pcol(:,j),'col');
% Next, PIE time propagator matrix
A_local=A_2PDEstate_local(1:Nrow-prow(i)+1,:);

 % Create a cell array for the full operator
         A_cell{i,j}=A_local;
         A_cell_2PDEstate{i,j}=A_2PDEstate_local;
     end
end

% Concatenating from cellular to a global matrix structure 

 A=cell2mat(A_cell);
 A_2PDEstate=cell2mat(A_cell_2PDEstate);