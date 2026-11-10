function zmap = Zmap(Zsrc,Ztgt)
% ZMAP = ZMAP(ZSRC,ZTGT) returns, for each monomial in the Kronecker
% product basis defined by the cell ZSRC of monomial degree arrays, the
% index of the same monomial in the Kronecker product basis defined by
% ZTGT.
%
% INPUTS
% - Zsrc:   1 x N cell, with element i a nmonsS(i) x deps(i) array, each
%           row specifying the degrees of one monomial along direction i;
% - Ztgt:   1 x N cell of the same form, with element i a nmonsT(i) x
%           deps(i) array. Each row of Zsrc{i} must appear in Ztgt{i};
%
% OUTPUTS
% - zmap:   nZS x 1 array of integers, for nZS = prod(nmonsS), with zmap(z)
%           the index in the target basis of monomial z of the source
%           basis;
%
% NOTES
% - In a Kronecker product basis, the first direction changes slowest.
%

% DJ, 10/15/2026: Initial coding;

N = numel(Zsrc);
[nS,stS,nZS] = Zdims(Zsrc);
[~,stT,~] = Zdims(Ztgt);
zmap = ones(nZS,1);
lin = (0:nZS-1)';
for i=1:N
    zi = mod(floor(lin/stS(i)),nS(i))+1;        % index of monomial in direction i
    [is_in,loc] = ismember(Zsrc{i},Ztgt{i},'rows');
    if ~all(is_in)
        error("Monomials of the source basis must be contained in the target basis.")
    end
    zmap = zmap + (loc(zi)-1)*stT(i);
end

end
