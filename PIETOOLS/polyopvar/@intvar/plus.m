function C = plus(A,B)
% C = PLUS(A,B) returns the 'intvar' object C representing the sum of the
% functionals defined by the 'intvar' objects A and B
%
% INPUTS
% - A:      m x n 'intvar' object acting on a degree-d distributed
%           monomial;
% - B:      m x n 'intvar' object acting on the same degree-d distributed 
%           monomial as A;
%
% OUTPUTS
% - C:      m x n 'intvar' object representing the sum of the functionals
%           defined by A and B;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus
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
% DJ, 02/03/2026: Initial coding

% Check that the operators are of the same dimension
Adim = A.matdim;    Bdim = B.matdim;
if any(Adim~=Bdim)
    error("Matrix dimensions of functionals must match for addition.")
end
% Check that the operators act on a distributed monomial of the same degree
% (we cannot check that they actually act on the same distributed monomial)
if size(A.omat,2)~=size(B.omat,2)
    error("Functionals must act on distributed monomial of same degree.")
end
% Check that the domain of the spatial variables matches
if any(any(A.dom~=B.dom))
    error("Domain of integration of the functionals must match.")
end

% Express the parameters defining the functional B in terms of the same 
% variables as those defining A
Avars = A.pvarname;     Aparams = A.params;
Bvars = B.pvarname;     Bparams = B.params;
if ~isequal(Avars,Bvars)
    p = size(Bvars,2);
    [Bvar_new,idx1,idx2] = unique([Bvars(1,:)'; Bparams.varname]);
    Bvar_new(idx1(idx1<=p)) = Avars(1,idx1(idx1<=p));        % replace Bvars by Avars
    Bparams.varname = Bvar_new(1,idx2(p+1:end));
end

% Concatenate the list of parameters and associated limits of integration
% of the functional B to that of A
params_full = [Aparams,Bparams];
omat_full = [A.omat; B.omat];

% Combine duplicate terms in the functional, corresponding to integrals
% with the same limits of integration
[P,omat] = uniquerows_integerTable(omat_full);
C = A;
C.omat = omat;
C.params = params_full*P;

end