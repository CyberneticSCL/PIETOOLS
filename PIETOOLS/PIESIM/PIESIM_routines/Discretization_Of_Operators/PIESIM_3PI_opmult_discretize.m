%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_opmult_discretize.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of a polymonial
% multiplicative operator component of 3PI operator
% Discretizes an operator component (NOT FULL OPERATOR)
%
%   (R v)(s) = Rop(s) v(s).
%
% Rop is sampled on the Chebyshev grid of degree deg+1 (deg = Rop.maxdeg)
% and transformed to coefficients a_j; column i+1 of A holds the
% coefficients of Rop T_i, up to degree N, from
% T_j T_i = (T_{i+j} + T_{|i-j|})/2, T_k the Chebyshev polynomials. Rop may
% depend on one variable, of any name; a polynomial in more variables
% throws an error.
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop -  polynomial multplicative operator R0 of size 1x1
%
% Outputs:
% A - square matrix of size N+1 x N+1 that represents a Chebyshev
% discretizaiton of the polynomial multiplicative matrix operator. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_3PI_opmult_discretize
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
% Initial coding YP  - 06_28_2022

function A=PIESIM_3PI_opmult_discretize(N, Rop)

A(1:N+1,1:N+1)=0;

% Degree of Rop, and its values on the Chebyshev grid of deg+2 nodes
if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(Rop,'polynomial') && ~isdouble(Rop)
    if numel(Rop.varname)>1
        error("The input polynomial depends on multiple variables; discretization not supported...")
    else
        var = polynomial(Rop.varname);
        Reval = subs(Rop,var,chebgrid);
    end    
else
    Reval=Rop*ones(1,deg+2);
end

% Chebyshev coefficients of Rop
acheb=fcht(double(Reval));

% Column i: coefficients of Rop*T_{i-1}, truncated at degree N
for i=1:N+1
    id=i-1;
    for j=1:length(acheb)
        jd=j-1;
        if (jd+id<=N)
        A(jd+id+1,i)=A(jd+id+1,i)+0.5*acheb(j);
        end
        if (abs(jd-id)<=N)
        A(abs(jd-id)+1,i)=A(abs(jd-id)+1,i)+0.5*acheb(j);
        end
    end
end




