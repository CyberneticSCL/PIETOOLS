%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_opint_discretize.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - a dicrete matrix representation of an integrative portion
% of the 3PI operator
% Discretizes an operator component (NOT FULL OPERATOR) 
%
%   (R v)(s) = int_{-1}^{s} R1(s,theta) v(theta) dtheta
%              + int_{s}^{1} R2(s,theta) v(theta) dtheta.
%
% R1 and R2 are sampled on the 2D Chebyshev grid of degree maxdeg+1 and
% transformed to coefficients R(s,theta) = sum_{j,m} c_{jm} T_j(s) T_m(theta),
% T_k the Chebyshev polynomials. Column k+1 of A holds the Chebyshev
% coefficients of R T_k, up to degree N, from
% T_m T_k = (T_{k+m} + T_{|k-m|})/2 and the antiderivatives of T_k.
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R1 - integrative operator R1 of size 1x1
% R2 - integrative operator R2 of size 1x1
%
% Output:
% A - a  matrix of size (N+1)x(N+1) that represents a Chebyshev
% discretizaiton of the integrative portion of the PI operator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_3PI_opint_discretize
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
% Initial coding YP  - 6_28_2022
function A=PIESIM_3PI_opint_discretize(N, R1, R2)

A(1:N+1,1:N+1)=0;

% Common degree of R1, R2, and their 2D Chebyshev coefficients
if isa(R1,'polynomial')
    deg1=R1.maxdeg;
else
    deg1=1;
end

if isa(R2,'polynomial')
    deg2=R2.maxdeg;
else
    deg2=1;
end

maxdeg=max(deg1,deg2);

Rcheb1=PIESIM_PI_discretize(R1,maxdeg);
Rcheb2=PIESIM_PI_discretize(R2,maxdeg);

for k=0:N
    % Finding a column vector of Chebyshev coefficient for the column k+1
    % of matrix A
vecmult1=zeros(N+1,1);
vecmult2=zeros(N+1,1);
for j=0:maxdeg+1
    for m=0:maxdeg+1
% Integrate R1 (T_(k+m) and T_|k-m|) polynomials on [-1,s]
% Integrate R2 (T_(k+m) and T_|k-m|) polynomials on [s,1]

vecint1=zeros(N+1,1);
vecint2=zeros(N+1,1);

i=[k+m,abs(k-m)];
for kk=i
 % int T_0 = T_1; the value at the fixed limit (-1 or 1) enters entry 1
 if (kk==0)
     vecint1(2)=vecint1(2)+0.5*Rcheb1(j+1,m+1);
     vecint2(2)=vecint2(2)-0.5*Rcheb2(j+1,m+1);
     vecint1(1)=vecint1(1)-0.5*Rcheb1(j+1,m+1)*(-1)^1;
     vecint2(1)=vecint2(1)+0.5*Rcheb2(j+1,m+1)*(1)^1;

  elseif (kk==1)
     % int T_1 = (T_2 + T_0)/4
     vecint1(1)=vecint1(1)+0.5*0.25*Rcheb1(j+1,m+1);
     vecint1(3)=vecint1(3)+0.5*0.25*Rcheb1(j+1,m+1);
     vecint2(1)=vecint2(1)-0.5*0.25*Rcheb2(j+1,m+1);
     vecint2(3)=vecint2(3)-0.5*0.25*Rcheb2(j+1,m+1);
     vecint1(1)=vecint1(1)-0.5*0.25*Rcheb1(j+1,m+1)*((-1)^0+(-1)^2);
     vecint2(1)=vecint2(1)+0.5*0.25*Rcheb2(j+1,m+1)*((1)^0+(1)^2);
else
    % int T_kk = (T_{kk+1}/(kk+1) - T_{kk-1}/(kk-1))/2, degrees above N
    % dropped
    if (kk+1+1<=N+1)
    vecint1(kk+1+1)=vecint1(kk+1+1)+0.5*0.5*Rcheb1(j+1,m+1)/(kk+1);
    vecint2(kk+1+1)=vecint2(kk+1+1)-0.5*0.5*Rcheb2(j+1,m+1)/(kk+1);
    end
    if (kk-1+1<=N+1)
    vecint1(kk-1+1)=vecint1(kk-1+1)-0.5*0.5*Rcheb1(j+1,m+1)/(kk-1);
    vecint2(kk-1+1)=vecint2(kk-1+1)+0.5*0.5*Rcheb2(j+1,m+1)/(kk-1);
    end
    vecint1(1)=vecint1(1)-0.5*0.5*Rcheb1(j+1,m+1)/(kk+1)*(-1)^(kk+1)+0.5*0.5*Rcheb1(j+1,m+1)/(kk-1)*(-1)^(kk-1);
    vecint2(1)=vecint2(1)+0.5*0.5*Rcheb2(j+1,m+1)/(kk+1)*(1)^(kk+1)-0.5*0.5*Rcheb2(j+1,m+1)/(kk-1)*(1)^(kk-1);
end % if
end % for kk=i

% Multiply by T_j(s)
   if (j==0) 
       vecmult1=vecmult1+vecint1;
       vecmult2=vecmult2+vecint2;
   else 
       % T_j T_i = (T_{i+j} + T_{|i-j|})/2, degrees above N dropped
       for i=1:N+1
           iorder=i-1;
       if (iorder+j+1<=N+1)
       vecmult1(iorder+j+1)=vecmult1(iorder+j+1)+0.5*vecint1(i);
       vecmult2(iorder+j+1)=vecmult2(iorder+j+1)+0.5*vecint2(i);
       end
       if (abs(iorder-j)+1<=N+1)
       vecmult1(abs(iorder-j)+1)=vecmult1(abs(iorder-j)+1)+0.5*vecint1(i);
       vecmult2(abs(iorder-j)+1)=vecmult2(abs(iorder-j)+1)+0.5*vecint2(i);
       end
       end % i
    end % j
    end % m
end % j
A(:,k+1)=vecmult1(:)+vecmult2(:);
end % k
