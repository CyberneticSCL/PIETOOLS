function [prog,Zop] = lpivar(prog,n,I,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Zop] = lpivar(prog,n,I,d) declares 
% a 4-PI operator variable of the form 
%
% P_{PQRS}[x] = [Px + \int_{I(1)}^(I(2))Q_1(s)y(s)ds
%         [y]   [Q_2(s)^T x + R_0(s)y(s)+\int_{I(1)}^(s) R_{1}(s,t)y(t)dt+\int_{s}^(I(2)) R_{2}(s,t)y(t)dt
% 
% INPUT 
%   prog: SOS program to modify.
%   n=[n11 n12; n21 n22]: dimension of the operator
%   I = [l u] spatial domain
%   -Optional INPUTS
%   d(1): degree of s in Z(s), translates to degree of var1 in Q1,Q2 and R0
%   d(2): degree of t in R1 and R2, defaults to d(1) if length of d=2,
%   d(3): degree of s in R1 and R2, defaults to d(2) if length of d=2,
%
% OUTPUT 
%   Zop: operator structure
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpivar
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019
%

if nargin<3
    error("Not enough input arguments");
elseif nargin==3
    d = [1,1,1];
end

pvar s theta;
var1 = s; var2 = theta;
Zs = monomials(var1,0:d(1));
[prog,P] = sospolymatrixvar(prog,monomials(var1,0),[n(1,1) n(1,2)]);
[prog,Q1] = sospolymatrixvar(prog,Zs,[n(1,1) n(2,2)]);
[prog,Q2] = sospolymatrixvar(prog,Zs,[n(2,1) n(1,2)]);
[prog,R0] = sospolymatrixvar(prog,Zs,[n(2,1) n(2,2)]);

if length(d)==2
    Zsth = mpmonomials({var1,var2},{0:d(2),0:d(1)});
elseif length(d)==1
    Zsth=mpmonomials({var1,var2},0:d);
elseif length(d)==3
    Zsth=mpmonomials({var1,var2},{0:d(3),0:d(2)});
else
    error('degree must be a scalar or a vector of length 2 or 3')
end
[prog,R1] = sospolymatrixvar(prog,Zsth,[n(2,1) n(2,2)]);
[prog,R2] = sospolymatrixvar(prog,Zsth,[n(2,1) n(2,2)]);

opvar Zop;
Zop.P = P;
Zop.Q1 = Q1; 
Zop.Q2 = Q2; 
Zop.R.R0 = R0; 
Zop.R.R1 = R1; 
Zop.R.R2 = R2; 
Zop.I = I; Zop.var1 = var1; Zop.var2 = var2;
Zop.dim = n;




