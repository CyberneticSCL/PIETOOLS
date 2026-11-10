function P = remap_cols(P,zmap,nZT)
% P = REMAP_COLS(P,ZMAP,NZT) expresses the sparse coefficient matrix P in
% a monomial basis of NZT monomials, by moving column (j-1)*nZS+z of P,
% associated with monomial z in matrix column j, to column
% (j-1)*NZT+ZMAP(z). A suitable map ZMAP is returned by 'Zmap' and
% 'Zpermmap'.
%
% INPUTS
% - P:      r x n*nZS sparse array, for nZS = numel(zmap), with column
%           (j-1)*nZS+z holding the coefficients of monomial z in matrix
%           column j;
% - zmap:   nZS x 1 array of integers in {1,...,nZT}, with zmap(z) the
%           index of monomial z in the new basis;
% - nZT:    scalar integer specifying the number of monomials in the new
%           basis;
%
% OUTPUTS
% - P:      r x n*nZT sparse array, with column (j-1)*nZT+z holding the
%           coefficients of monomial z of the new basis in matrix column j;
%
% NOTES
% - If ZMAP maps several monomials to the same monomial, their
%   coefficients are added.
%

% DJ, 10/15/2026: Initial coding;

nZS = numel(zmap);
if nZT==nZS && isequal(zmap,(1:nZS)')
    return
end
r = size(P,1);
n = size(P,2)/nZS;
if n==1 && isequal(zmap,(1:nZS)')
    % Only append empty columns
    P = [P,sparse(r,nZT-nZS)];
    return
end
[i,j,v] = find(P);
i = i(:);   j = j(:);   v = v(:);       % row vectors if P has one row
jb = ceil(j/nZS);
z = j - (jb-1)*nZS;
P = sparse(i,(jb-1)*nZT+zmap(z),v,r,n*nZT);

end
