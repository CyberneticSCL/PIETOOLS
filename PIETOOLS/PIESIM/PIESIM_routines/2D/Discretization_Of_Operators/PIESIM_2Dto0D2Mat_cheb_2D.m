%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_2Dto0D2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of a 2D to 0D
% operator (R02)
%
%   (R v) = int_{-1}^{1} int_{-1}^{1} R(x,y) v(x,y) dx dy.
%
% A monomial term c*x^a*y^b of R(i,j) contributes kron(I_y, I_x), with I_x
% the row of integrals int_{-1}^{1} c*x^a*T_k(x) dx, k = 0,...,N_x-p(1,j),
% and I_y that of y^b, k = 0,...,N_y-p(2,j); T_k the Chebyshev polynomials.
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% R -  polymonial block corresponding to R20 operator
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% p - 2 x ns matrix with the degrees of differentiability in x and y for
%     the 2D states
%
%
% Outputs:
% A - Chebyshev discretizaiton of R02 block that represents a block of a
% total matrix operator for time-advancement of the
% spatially-discretized PIE solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_2Dto0D2Mat_cheb_2D
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


function A=PIESIM_2Dto0D2Mat_cheb_2D(N, R, var1, p)

% Extract second spatial variable (y) defining the operator R
s2 = var1(2);

% Numbers of 2D states and of outputs
ns=size(p,2);
no=size(R,1);

% Row i of A: the blocks of all 2D states j, each the sum over the
% monomial terms of R(i,j)
for i=1:no
    Astring=zeros(1,0);
    for j=1:ns
        Rloc=R(i,j);
        int=0;
    for k=1:Rloc.nterms
        % Term k of R(i,j); y = 1 leaves its factor c*x^a
        Rstrip=polynomial(Rloc.coeff(k),Rloc.degmat(k,:),Rloc.varname,Rloc.matdim);
            Rs1=subs(Rstrip,s2,1);
            % Integrate in s1 direction over the interval [-1,1]
            int_s1=PIESIM_PI2Mat_opint_cheb(N(1), Rs1, p(1,j));
            % Factor y^b of term k
            index = find(strcmp(Rstrip.varname, s2.varname{1}));
            Rs2=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            % Integrate in s2 direction over the interval [-1,1]
            int_s2=PIESIM_PI2Mat_opint_cheb(N(2), Rs2, p(2,j));
            % Add the Kronecker product of the two rows (x index fastest)
            int=int+kron(int_s2,int_s1);
    end
            % Append block j
            Astring=[Astring,int];
    end
    A(i,:)=Astring;
end




