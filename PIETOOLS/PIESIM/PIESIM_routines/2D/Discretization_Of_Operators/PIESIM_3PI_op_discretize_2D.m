%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_op_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of 1D to 2D
% operators (R2x, R2y)
%
%   flag = 0:      (R v)(x,y) = R(x,y) v(x),
%   flag = [1 0]:  (R v)(x,y) = int_{-1}^{x} R(x,y,theta) v(theta) dtheta,
%   flag = [0 1]:  (R v)(x,y) = int_{x}^{1} R(x,y,theta) v(theta) dtheta,
%
% for R2x (dir = 'x'), and likewise with x and y exchanged for R2y. A
% monomial term of R contributes the Kronecker product of the Chebyshev
% coefficients of its factor in y with the square 1D discretization of its
% factor in (x,theta), x index fastest. A has (N_x+1)(N_y+1) rows and
% N_s+1 columns (s = dir).
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% R -  polymonial block (one of the components of the 3PI operator -
% multiplicative or integraitve) of R2x or R2y
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% dir - direction ('x' for R2x and 'y' for R2y)
% flag - determines if R0, R1, or R2 component of 3PI operator is given
% flag = 0 - multiplicative (R0) operator
% flag = [1 0] - integrative R1 operator
% flag = [0 1] - integrative R2 operator
%
%
% Outputs:
% A - Chebyshev
% discretizaiton of R0, R1 or R2 component of the R2x/R2y operator
% (depending on the flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_3PI_op_discretize_2D
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

function A=PIESIM_3PI_op_discretize_2D(N, R, var1, dir, flag)

% Extract spatial variables (x,y) defining the operator R
s1 = var1(1);   s2 = var1(2);

% Variable other than dir (var); orders Ncol in dir and Nrow in var
if dir=='x'
    var=s2;
    snstr=s2.varname{1};
    Ncol=N(1);
    Nrow=N(2);
else
    var=s1;
    snstr=s1.varname{1};
    Ncol=N(2);
    Nrow=N(1);
end

% Sum over the monomial terms of R
A=zeros(prod(N+1),Ncol+1);

    for k=1:R.nterms
        % Term k of R; var = 1 leaves its factor in dir
        Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
            Rcomp=subs(Rstrip,var,1);
            if (flag==0)
            % Multiplicative operator in "dir" direction 
            op_comp=PIESIM_3PI_opmult_discretize(Ncol, Rcomp);
            else 
            % Integrative 3PI operator in "dir" direction 
            Rcomp=Rcomp*flag;
            op_comp=PIESIM_3PI_opint_discretize(Ncol, Rcomp(1), Rcomp(2));
            end
            % Poly2mat in the other direction
            index = find(strcmp(Rstrip.varname, snstr));
            Rpoly=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            op_poly=PIESIM_Poly2Mat_cheb(Nrow, Rpoly, 0);
            % Add the Kronecker product column by column (x index fastest)
            if (dir=='x')
            for jcol=1:Ncol+1
            A(:,jcol)=A(:,jcol)+kron(op_poly,op_comp(:,jcol));
            end
            else
            for jcol=1:Ncol+1
            A(:,jcol)=A(:,jcol)+kron(op_comp(:,jcol),op_poly);
            end
            
            end
    end % k