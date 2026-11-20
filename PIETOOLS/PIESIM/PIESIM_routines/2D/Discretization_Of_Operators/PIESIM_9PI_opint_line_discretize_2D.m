%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI_opint_line_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of
% singly integrated blocks (R01, R02, R10, R20) of 9PI operator in 2D
% Discretizes an operator component of 9PI operator (NOT FULL OPERATOR)
%
%   sign = -1:  (R v)(x,y) = int_{-1}^{x} R(x,y,theta) v(theta,y) dtheta,
%   sign =  1:  (R v)(x,y) = int_{x}^{1} R(x,y,theta) v(theta,y) dtheta,
%
% for dir = 'x', and likewise in y for dir = 'y'. Each monomial term of R
% splits into a factor in (x,theta), discretized as a 1D integral
% operator, and a factor in y, discretized as a 2D multiplier. A maps the
% 2D Chebyshev coefficients of v to those of R v, x index fastest.
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% R -  polymonial multiplicative block (R01, R02, R10 or R20) of 9PI operator in 2D
% corresponding to a single solution state
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% var2: 2x1 pvar array specifying the dummy variables for integration
%           (theta,nu)
% dir - direction of integration ('x' or 'y')
% sign  - (-1 for [-1,x] or [-1,y] integration and 1 for [x, 1] or [y, 1] integation) 
%
% Outputs:
% A - Chebyshev
% discretizaiton of R01, R02, R10 or R20 block that represents a block of a
% total matrix operator for time-advancement of the
% spatially-discretized PIE solution (spatially-discrertized propagator
% matrix) - square prod(N+1) x prod(N+1) matrix, not truncated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_9PI_opint_line_discretize_2D
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

function A=PIESIM_9PI_opint_line_discretize_2D(N, R, var1, var2, dir, sign)

% Extract spatial and dummy variables defining the operator R
s1 = var1(1);   s1_dum = var2(1);
s2 = var1(2);   s2_dum = var2(2);

% snative, snstr correspond to the variable not affected by integration (not in the differential or limits) 
% s is variable in limits, theta is the dummy variable of integration (differential)
R_native=false;

if (dir=='x')
    s=s1;
    theta=s1_dum;
    snative=s2;
    snstr=s2.varname{1};
    idir=1;
else
    s=s2;
    theta=s2_dum;
    snative=s1;
    snstr=s1.varname{1};
    idir=2;
end

% Degree of R, and the Chebyshev grid of deg+2 nodes
if isa(R,'polynomial')
    deg=R.maxdeg;
else
    deg=1;
end

% Values of R on the grid in (theta,s): one array per monomial term if R
% depends on snative (snative = 1 in each), else one array
chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(R,'polynomial')
        check_native=ismember(R.varname,snstr);
        if (any(check_native))
            R_native=true;
            for k=1:R.nterms
            Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
            Rstrip=subs(Rstrip,snative,1);
            Reval{k}=subs(subs(Rstrip,theta,chebgrid),s,chebgrid); 
            end
        else
        Reval{1}=subs(subs(R,theta,chebgrid),s,chebgrid); 
        end
    else % not a polymoial 
    Reval{1}=R*ones(deg+2);
    end


% 2D Chebyshev coefficients of each array
for k=1:numel(Reval)
Revald{k}=double(Reval{k});
Rchebcell{k}=fcgltran2d(Revald{k},1);
end

% Sum over the arrays of the 1D integral operator in dir, composed with the
% multiplier by the factor in snative
for kcell=1:numel(Reval)
    Rcheb=Rchebcell{kcell};

    Aint=PIESIM_integral_projection1D(N(idir), Rcheb, sign);

% Construct multiplicative operator in the other direction

    Multop=eye(prod(N+1));
    Idx=eye(N(2)+1);
    Idy=eye(N(1)+1);
   
if (R_native)
    index = find(strcmp(R.varname, snstr));
    Rop=polynomial(1,R.degmat(kcell,index),R.varname(index),R.matdim);
    Multop=PIESIM_9PI_opmult_discretize_2D(N, Rop, var1);
end

% Identity in the direction other than dir (x index fastest)
if (dir=='x')
    if (exist('A','var'))
    A=A+kron(Idx,Aint)*Multop;
   else
    A=kron(Idx,Aint)*Multop;
    end
else
      if (exist('A','var'))
    A=A+kron(Aint,Idy)*Multop;
     else
    A=kron(Aint,Idy)*Multop;
      end
end

end % kcell




