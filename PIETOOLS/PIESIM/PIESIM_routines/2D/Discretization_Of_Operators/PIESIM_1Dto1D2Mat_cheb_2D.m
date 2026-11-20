%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_1Dto1D2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_2PDEstate - discrete matrix representations of 1D to 1D
% cross-state operators (Rxy, Ryx)
% Discretized FULL OPERATOR of 1D-to-1D states
%
%   (R v)(t) = int_{-1}^{1} R(s,t) v(s) ds,
%
% with s the direction of integration and t the other direction. A
% monomial term c*s^a*t^b of R(i,j) contributes the outer product of the
% Chebyshev coefficients of t^b with the row of integrals
% int_{-1}^{1} c*s^a*T_k(s) ds, k = 0,...,N_s, T_k the Chebyshev polynomials.
% Block (i,j) of A_2PDEstate keeps the first N_s-pcol(j)+1 columns; block
% (i,j) of A keeps, of these, the first N_t-prow(i)+1 rows.
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% R -  polymonial block corresponding to Rxy or Ryx operator
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% prow - vector with the degrees of differentiability of the 1D row
%        (output) states, which depend on t only
% pcol - vector with the degrees of differentiability of the 1D column
%        (input) states, which depend on s only
% dir - direction of integration ('y' for Rxy and 'x' for Ryx)
%
%
% Outputs:
% A - Chebyshev
% discretizaiton of Rxy or Ryx block that represents a block of a
% total matrix operator for time-advancement of the
% spatially-discretized PIE solution 
% A_2PDEstate - Chebyshev
% discretizaiton of Rxy or Ryx block that represents a block of a
%  total matrix operator for reconstruction of the primary (PDE)
% solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_1Dto1D2Mat_cheb_2D
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

function [A, A_2PDEstate]=PIESIM_1Dto1D2Mat_cheb_2D(N, R, var1, prow, pcol, dir)

% Extract spatial variables (x,y) defining the operator R
s1 = var1(1); s2 = var1(2);

% Variable t (var) other than the direction of integration; orders Nrow in
% t and Ncol in s
if dir=='x'
    var=s2;
    snstr=s2.varname{1};
    Nrow=N(2);
    Ncol=N(1);
else
    var=s1;
    snstr=s1.varname{1};
    Nrow=N(1);
    Ncol=N(2);
end

% Numbers of row and column states
    ns_row=size(prow,2);
    ns_col=size(pcol,2);

% Discretize each block R(i,j) as the sum over its monomial terms
for i=1:ns_row
    rsize=Nrow-prow(i)+1;
    for j=1:ns_col
        Rloc=R(i,j);
        csize=Ncol-pcol(j)+1;
        int=zeros(Nrow+1,Ncol+1);
    for k=1:Rloc.nterms
        % Term k of R(i,j); t = 1 leaves its factor c*s^a
        Rstrip=polynomial(Rloc.coeff(k),Rloc.degmat(k,:),Rloc.varname,Rloc.matdim);
            Rint=subs(Rstrip,var,1);
            % Integrate in pint direction over the interval [-1,1]
            op_int=PIESIM_PI2Mat_opint_cheb(Ncol, Rint, 0);
            % Factor t^b of term k
            index = find(strcmp(Rstrip.varname, snstr));
            Rrow=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            % Convert to a matrix operator in prow direction 
            op_row=PIESIM_Poly2Mat_cheb(Nrow, Rrow, 0);
            % Add the outer product of the two factors
            int=int+kron(op_row,op_int);
    end % k

    % Truncate the total operator to the correct dimensions

% First, PIE-to-PDE state matrix
     int_2PDEstate=int(:,1:csize);
% Next, PIE time propagator matrix
     int_propagator=int_2PDEstate(1:rsize,:);

    % Store block (i,j)
    A_cell{i,j}=int_propagator;
    A_cell_2PDEstate{i,j}=int_2PDEstate;
    end % j
end %i

% Assemble the blocks into the global matrices
 A=cell2mat(A_cell);
 A_2PDEstate=cell2mat(A_cell_2PDEstate);







