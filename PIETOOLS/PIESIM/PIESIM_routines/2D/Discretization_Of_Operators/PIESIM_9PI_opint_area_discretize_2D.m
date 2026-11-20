%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI_opint_area_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of
% double-integrated blocks (R11, R12, R21, R22) of 9PI operator in 2D
% Discretizes operator component of 9PI operator (NOT FULL OPERATOR)
%
%   (R v)(x,y) = int_{lim(1,:)} int_{lim(2,:)} R(x,y,theta,nu) v(theta,nu) dnu dtheta,
%
% with lim(1,:) = [-1 x] or [x 1] and lim(2,:) = [-1 y] or [y 1]. R is
% sampled on the 4D Chebyshev grid of degree deg+1 (deg = R.maxdeg) and
% transformed to Chebyshev coefficients; the integral in x is projected for
% each pair of y-coefficients, then the integral in y. A maps the 2D
% Chebyshev coefficients of v to those of R v, x index fastest.
%
% Inputs:
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% R -  polymonial block (R11, R12, R21 or R22) of 9PI operator in 2D
% corresponding to a single solution state
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% var2: 2x1 pvar array specifying the dummy variables for integration
%           (theta,nu)
% lim  - limits of 2D integraton (for example, [-1 sx;-1 sy]) 
%
% Outputs:
% A - Chebyshev
% discretizaiton of R11, R12, R21 or R22 block that represents a block of a
% total matrix operator - size prod(N+1) x prod(N+1), not truncated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_9PI_opint_area_discretize_2D
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

function A=PIESIM_9PI_opint_area_discretize_2D(N, R, var1, var2, lim)

% Extract spatial and dummy variables defining the operator R
s1 = var1(1);   s1_dum = var2(1);
s2 = var1(2);   s2_dum = var2(2);

% Degree of R, and the Chebyshev grid of deg+2 nodes
if isa(R,'polynomial')
    deg=R.maxdeg;
else
    deg=1;
end

Reval=zeros(deg+2,deg+2,deg+2,deg+2);

chebgrid=cos(pi*(0:deg+1)/(deg+1));

% Values of R on the 4D grid: substitute nu, then y, then theta and x,
% replicating along each variable R does not depend on
if isa(R,'polynomial')
% Rynu(i,j): R at y = chebgrid(i), nu = chebgrid(j)
if ismember(s2_dum.varname{1},R.varname)
    Ry=subs(R,s2_dum,chebgrid);
else
    Ry=R*ones(1,deg+2);
end
if ismember(s2.varname{1},R.varname)
    Rynu=subs(Ry,s2,chebgrid);
else
    Rynu=repmat(Ry,deg+2,1);
end
% Reval(:,:,i,j): x along the rows, theta along the columns
if ismember(s1_dum.varname{1},R.varname) && ismember(s1.varname{1},R.varname)
    for j=1:deg+2
        for i=1:deg+2
        Rmat=subs(subs(Rynu(i,j),s1_dum,chebgrid),s1,chebgrid);
        Reval(:,:,i,j)=Rmat(:,:);
        end
    end
elseif ismember(s1.varname{1},R.varname)
   % R independent of theta: replicated along the columns
   for j=1:deg+2
        for i=1:deg+2
        Re=subs(Rynu(i,j),s1,chebgrid);
        Reval(:,:,i,j)=repmat(Re',1,deg+2);
        end
   end
   else
       % R independent of x: replicated along the rows
       for j=1:deg+2
        for i=1:deg+2
       Re=subs(Rynu(i,j),s1_dum,chebgrid);
       Reval(:,:,i,j)=repmat(Re,deg+2,1);
        end
       end
   end
else
    % Constant R
    Reval=R*ones(deg+2,deg+2,deg+2,deg+2);
end

% 4D Chebyshev coefficients of R
Reval=double(Reval);
Rcheb=fcgltran4d(Reval);

% Integral in x over lim(1,:), for each pair (p,q) of y-coefficients
for q=1:deg+2
    for p=1:deg+2
        A(1:N(1)+1,1:N(1)+1,p,q)=PIESIM_integral_projection1D(N(1), Rcheb(:,:,p,q), lim(1,:));
    end
end

% Integral in y over lim(2,:), for each entry (i,k) of the x-projection
 for k=1:N(1)+1
     for i=1:N(1)+1
         Rycheb(:,:)=A(i,k,:,:);
         A2D(i,k,1:N(2)+1,1:N(2)+1)=PIESIM_integral_projection1D(N(2), Rycheb, lim(2,:));
     end
 end

% Order the rows and columns as 2D coefficient vectors, x index fastest
acheb=permute(A2D,[1 3 2 4]);
A=reshape(acheb,prod(N+1),prod(N+1));




