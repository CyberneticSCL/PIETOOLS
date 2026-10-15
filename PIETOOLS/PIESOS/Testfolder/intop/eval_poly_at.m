function K = eval_poly_at(P,vals)
% K = EVAL_POLY_AT(P,VALS) evaluates the 'polynomial' (or 'double') matrix
% P at the numeric values VALS of all of its variables, returning the
% resulting numeric matrix.
%
% INPUTS
% - P:      m x n 'polynomial' or 'double' object;
% - vals:   containers.Map from variable name (as in P.varname) to a
%           numeric value;
%
% OUTPUTS
% - K:      m x n 'double' array, the value of P at VALS. If P is a
%           'double', it is returned unchanged;
%
% NOTES
% - Monomials are the rows of P.degmat, and the coefficients P.C have size
%   nZ x m*n, with column (j-1)*m+i corresponding to entry (i,j), as in
%   eval_dpvar_at.m.

% DJ (CC), 10/15/2026: Initial coding;

if isa(P,'double')
    K = P;
    return
end
[m,n] = size(P);
degmat = full(P.degmat);
nZ = size(degmat,1);
z = ones(nZ,1);
for v=1:numel(P.varname)
    z = z.*vals(P.varname{v}).^degmat(:,v);
end
K = reshape(full(z.'*P.C),m,n);

end
