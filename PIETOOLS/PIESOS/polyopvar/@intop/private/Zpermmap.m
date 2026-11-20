function zmap = Zpermmap(Zin,perm)
% ZMAP = ZPERMMAP(ZIN,PERM) returns, for each monomial in the Kronecker
% product basis defined by the cell ZIN of arrays of monomial degrees, the
% index of the same monomial in the Kronecker product basis defined by
% ZIN(PERM), in which the order of the variables is changed, so that
% variable PERM(p) becomes variable p.
%
% INPUTS
% - Zin:    1 x N cell, with element i a nmons(i) x deps(i) array, each row
%           specifying the degrees of one monomial along direction i;
% - perm:   1 x N array of integers, a permutation of 1,...,N;
%
% OUTPUTS
% - zmap:   nZ x 1 array of integers, for nZ = prod(nmons), with zmap(z)
%           the index in the basis defined by ZIN(PERM) of monomial z of
%           the basis defined by ZIN;
%
% NOTES
% - The map is suitable as input to 'remap_cols'.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Zpermmap
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 10/15/2026: Initial coding;

N = numel(Zin);
[nS,~,nZ] = Zdims(Zin);
Ztgt = Zin(perm);
[~,stT,~] = Zdims(Ztgt);
stS = zeros(1,N);
for i=1:N
    stS(i) = prod(nS(i+1:end));
end
lin = (0:nZ-1)';
zmap = ones(nZ,1);
for p=1:N
    i = perm(p);
    zi = mod(floor(lin/stS(i)),nS(i))+1;        % index of the monomial along source variable i
    zmap = zmap + (zi-1)*stT(p);
end

end
