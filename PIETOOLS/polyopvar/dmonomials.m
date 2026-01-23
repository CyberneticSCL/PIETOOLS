function Z = dmonomials(vartab,degs)
% Z = MONOMIALS(VARTAB,DEGS) constructs a vector of distributed monomials 
% of degrees DEGS in the variables VARTAB.
%
% INPUTS
% - vartab: p x 1 array of type 'polyopvar', specifying the variables in
%           which to construct a vector of monomials;
% - degs:   n x 1 array specifying the cumulative degrees of the monomials
%           to be included in the vector;
%
% OUTPUTS
% - Z:      m x 1 array of type 'polyopvar' representing a vector of all
%           monomials of cumulative degrees specified by "degs";
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dmonomials
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

if ~isa(vartab,'polyopvar')
    error("Variables to appear in distributed monomials must be specified as object of type 'polyopvar'.")
end

% Check that the degrees make sense
if any(degs<0) || any(round(degs)~=degs)
    error("Monomial degrees must be specified as nonnegative integers");
end

% Build the matrix of cumulative degrees at most degs in the variables
nvars = numel(vartab.varname);
degmat = zeros(0,nvars);
for ii=1:numel(degs)
    d = degs(ii);
    degmat_ii = sparse(1,nvars);
    for jj = 1:nvars
        nZ = size(degmat_ii,1);
        degmat_ii = sprepmat(degmat_ii,d+1,1);
        for kk = 0:d
            degmat_ii(nZ*kk+(1:nZ),jj) = kk;
        end
        degmat_ii = degmat_ii(sum(degmat_ii,2)<=d,:);   % Throw away invalid monomials
    end
    degmat = [degmat;degmat_ii(sum(degmat_ii,2)==d,:)];
end

% Declare the distributed monomial basis
Z = vartab;
Z.degmat = degmat;       

end