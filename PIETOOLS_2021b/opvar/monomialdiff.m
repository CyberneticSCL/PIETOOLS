function [d1, d2] = monomialdiff(Z1,Z2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [d1, d2] = monomialdiff(Z1,Z2) function takes two polynomials and checks if they have same monomials.
% 
% INPUT
% Z1, Z2 : polynomials
% 
% OUTPUT
% d1: monomials in Z1 that are not in Z2
% d2: monomials in Z2 that are not in Z1
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PIETools - monomialdiff
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS  - 7_26_2019
%

if ~isa(Z1,'polynomial')||~isa(Z2,'polynomial')
    error('Inputs must be polynomials');
end

var1 = Z1.varname; var2 = Z2.varname; 
[diff1, idx1] = setdiff(var1,var2);
[diff2, idx2] = setdiff(var2,var1);

a = full(Z1.degmat);
b = full(Z2.degmat);

a(:,end+1:end+length(idx2)) = 0; 
b(:,end+1:end+length(idx1)) = 0; 

d1.varname = [var1, diff2]; d2.varname = [var2, diff1];
[~, idx1] = sort(d1.varname);[sortedvar, idx2] = sort(d2.varname);
degmat1 = a(:,idx1); degmat2 = b(:,idx2);

[L1,~] = ismember(degmat1,degmat2,'rows');  [L2,~] = ismember(degmat2,degmat1,'rows');   
d1.degmat = degmat1(~L1,:); d1.coefficient = speye(size(d1.degmat,1)); d1.matdim = [size(d1.degmat,1),1];
d2.degmat = degmat2(~L2,:); d2.coefficient = speye(size(d2.degmat,1)); d2.matdim = [size(d2.degmat,1),1];

d1 = polynomial(d1.coefficient,d1.degmat,sortedvar,d1.matdim);
d2 = polynomial(d2.coefficient,d2.degmat,sortedvar,d2.matdim);
end