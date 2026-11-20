%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI_opmult_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of a polymonial
% multiplicative block (R00) of 9PI operator in 2D
% Discretizes an operator component (not a full operator)
%
%   (R v)(x,y) = R(x,y) v(x,y).
%
% A monomial term c*x^a*y^b of R contributes kron(M_y, M_x), with M_x, M_y
% the 1D Chebyshev multiplication matrices of c*x^a and y^b (x index
% fastest).
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% R -  polymonial multiplicative block (R00) of 9PI operator in 2D
% corresponding to a single solution state
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
%
% Outputs:
% A - Chebyshev
% discretizaiton of a R00 block that represents a block of a
% spatially-discretized total matrix operator - size prod(N+1) x prod(N+1),
% not truncated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_9PI_opmult_discretize_2D
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

function A=PIESIM_9PI_opmult_discretize_2D(N, R, var1)

% Extract second spatial variable (y) defining the operator R
s2 = var1(2);

% Sum over the monomial terms of R
A=0;
for k=1:R.nterms
    % Term k of R; y = 1 leaves its factor c*x^a
    Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
        Rs1=subs(Rstrip,s2,1);
        % Multiplicative operator in s1 direction
        mult_s1=PIESIM_3PI_opmult_discretize(N(1), Rs1);
        % Factor y^b of term k
        index = find(strcmp(Rstrip.varname, s2.varname{1}));
        Rs2=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
        % Multiplicative operator in s2 direction
        mult_s2=PIESIM_3PI_opmult_discretize(N(2), Rs2);
        % Add the Kronecker product (x index fastest)
        A=A+kron(mult_s2,mult_s1);
end