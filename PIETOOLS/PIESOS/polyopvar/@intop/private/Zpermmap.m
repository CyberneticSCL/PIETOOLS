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
