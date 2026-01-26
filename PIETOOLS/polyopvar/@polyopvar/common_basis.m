function [A_out,B_out] = common_basis(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_BASIS(A_IN,B_IN) takes two distributed polynomials
% and expresses them in terms of the same basis of distributed monomials
%
% INPUTS
% - A_in, B_in:     m x n 'polyopvar' objects
%
% OUTPUTS
% - A_out, B_out:   m x n 'polyopvar' objects representing the same
%                   distributed polynomials as A_in, B_in, respectively,
%                   but now expressed in terms of the same basis of
%                   monomials, so that
%                   A_out.varname = B_out.varname;
%                   A_out.degmat = B_out.degmat;
%                   A_out.pvarname = B_out.pvarname;
%                   A_out.varmat = B_out.varmat;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_basis
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
% DJ, 01/22/2026: Initial coding

% First, make sure the operators are expressed in terms of the same state
% variables and independent variables
[A,B] = common_vars(A_in,B_in);
degs1 = A.degmat;       degs2 = B.degmat;
C_A = A.C;              C_B = B.C;

% Combine the monomial bases into a unique basis
nZ1 = size(degs1,1);
degs_full = [degs1;degs2];
[Pmat,degs_new] = uniquerows_integerTable(degs_full);   % Pmat*Znew = Z_full
nZ = size(degs_new,1);
new2old_Cdcs1 = Pmat(1:nZ1,:)*(1:nZ)';
new2old_Cdcs2 = Pmat(nZ1+1:end,:)*(1:nZ)';
C1 = tensopvar(1,nZ);      C2 = tensopvar(1,nZ);
C1(1,new2old_Cdcs1) = C_A;
C2(1,new2old_Cdcs2) = C_B;

% Build polynomials in terms of the shared basis
A_out = A;                      B_out = B;
A_out.C = C1;                   B_out.C = C2;
A_out.degmat = degs_new;        B_out.degmat = degs_new;

end