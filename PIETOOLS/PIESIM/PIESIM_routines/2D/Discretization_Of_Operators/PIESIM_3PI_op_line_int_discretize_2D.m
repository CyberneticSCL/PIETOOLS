 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_op_line_int_discretize_2D.m    PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of single-component operators
% (line integral of 3PI operators) entering 2D to 1D operators (Rx2, Ry2)
% Dicretizes R0, R1 or R2 components (NOT FULL OPERATOR) depending on the
% `flag' variable
%
%   flag = 0:     (R v)(x) = int_{-1}^{1} R(x,y) v(x,y) dy,
%   flag = [1 0]: (R v)(x) = int_{-1}^{1} int_{-1}^{x} R(x,theta,y) v(theta,y) dtheta dy,
%   flag = [0 1]: (R v)(x) = int_{-1}^{1} int_{x}^{1} R(x,theta,y) v(theta,y) dtheta dy,
%
% for Rx2 (dir = 'x'), and likewise with x and y exchanged for Ry2. A
% monomial term of R contributes the Kronecker product of the row of
% integrals int_{-1}^{1} y^b T_k(y) dy, k = 0,...,N_y, with each row of the
% square 1D discretization of its factor in (x,theta), x index fastest;
% T_k the Chebyshev polynomials. A has N_s+1 rows (s = dir) and
% (N_x+1)(N_y+1) columns.
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% R -  polymonial block (one of the components of the 3PI operator -
% multiplicative or integraitve) of Rx2 or Ry2
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% dir - direction ('x' for Rx2 and 'y' for Ry2)
% flag - determines if R0, R1, or R2 component of 3PI operator is given
% flag = 0 - multiplicative (R0) operator
% flag = [1 0] - integrative R1 operator
% flag = [0 1] - integrative R2 operator
%
%
% Outputs:
% A - Chebyshev discretizaiton of the R0, R1 or R2 component of the 3PI
% operator integrated in the `dir' direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_3PI_op_line_int_discretize_2D
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

function A=PIESIM_3PI_op_line_int_discretize_2D(N, R, var1,dir, flag)

% Extract spatial variables (x,y) defining the operator R
s1 = var1(1);   s2 = var1(2);

% Variable of integration other than dir (var); orders Nrow in dir and
% Ncol in var
if dir=='x'
    var=s2;
    snstr=s2.varname{1};
    Nrow=N(1);
    Ncol=N(2);
else
    var=s1;
    snstr=s1.varname{1};
    Nrow=N(2);
    Ncol=N(1);
end
% Sum over the monomial terms of R
        A=zeros(Nrow+1,prod(N+1));

    for k=1:R.nterms
        % Term k of R; var = 1 leaves its factor in dir
        Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
            Rcomp=subs(Rstrip,var,1);
            if (flag==0)
            % Multiplicative operator in "dir" direction 
            op_comp=PIESIM_3PI_opmult_discretize(Nrow, Rcomp);
            else 
            % Integrative 3PI operator in "dir" direction 
            Rcomp=Rcomp*flag;
            op_comp=PIESIM_3PI_opint_discretize(Nrow, Rcomp(1), Rcomp(2));
            end
            % Integrate in the other direction over the interval [-1,1]
            index = find(strcmp(Rstrip.varname, snstr));
            Rint=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            op_int=PIESIM_PI2Mat_opint_cheb(Ncol, Rint, 0);
            % Add the Kronecker product row by row (x index fastest)
            if (dir=='x')
            for jrow=1:Nrow+1
            A(jrow,:)=A(jrow,:)+kron(op_int,op_comp(jrow,:));
            end
            else
            for jrow=1:Nrow+1
            A(jrow,:)=A(jrow,:)+kron(op_comp(jrow,:),op_int);
            end
            end
    end % k



