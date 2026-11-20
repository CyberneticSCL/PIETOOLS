%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI_discretize.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs Rcheb - PI operator evaluated on a Chebyshev grid
%
% The kernel R(s,theta) is evaluated on the (maxdeg+2) x (maxdeg+2) grid of
% Chebyshev nodes cos(pi*k/(maxdeg+1)), k = 0,...,maxdeg+1, and transformed
% to the coefficients Rcheb(j+1,m+1) of T_j(s) T_m(theta), T_k the
% Chebyshev polynomials. A kernel in one variable is taken as a function
% of s if the name of that variable has at most 2 characters, else of
% theta, and is replicated along the other.
%
% Inputs:
% R - compnent PI operator to be discretized
% maxdeg - maximum polynomial degree of a PI operator R
%
% Output:
% Rcheb - a square matrix of size (maxdeg+2)x(maxdeg+2) that evaluates the PI operator
% on a Chebysheb grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_PI_discretize
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
% Initial coding YP  - 12_21_2024
function Rcheb=PIESIM_PI_discretize(R,maxdeg)

% Values of R on the Chebyshev grid
chebgrid=cos(pi*(0:maxdeg+1)/(maxdeg+1));
if isa(R,'polynomial')
    switch(R.nvars)
        case 2
        % Two variables present
            Reval=subs(subs(R,R.varname(2),chebgrid),R.varname(1),chebgrid);
        case 1
        % One variable present
            Re=subs(R,R.varname,chebgrid);
            var=cell2mat(R.varname);
            if (length(var)<=2)
                Reval=repmat(Re',1,maxdeg+2);
            else
                Reval=repmat(Re,maxdeg+2,1);
            end
        case 0 
            % Scalar polynomial
            Reval=R*ones(maxdeg+2);
end
else
    Reval=R*ones(maxdeg+2);
end

% Nodal-to-spectral 2D Chebyshev transform
Reval=double(Reval);
Rcheb=fcgltran2d(Reval,1);