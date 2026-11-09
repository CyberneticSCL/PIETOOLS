function Ppoly = eval_dpvar_at(P,d0map)
% PPOLY = EVAL_DPVAR_AT(P,D0MAP) evaluates the 'dpvar' object P at fixed
% numeric values for its decision variables, returning the resulting
% 'polynomial'. Used as ground truth for testing 'intop' objects with
% decision-variable kernels.
%
% INPUTS
% - P:      m x n 'dpvar', 'polynomial', or 'double' object;
% - d0map:  containers.Map from decision-variable name (as in P.dvarname)
%           to a numeric value. Values are looked up by name, so P may
%           carry a reordered or partial subset of a larger, shared
%           decision-variable list;
%
% OUTPUTS
% - Ppoly:  m x n 'polynomial' object, equal to P with its decision
%           variables replaced by the values in d0map. If P is not a
%           'dpvar', it is returned unchanged;
%
% NOTES
% - dpvar layout: rows are matrix row outer, decision variable inner
%   (nd+1 rows per matrix row, constant term first); columns are matrix
%   column outer, monomial inner. 'polynomial' uses monomial rows and
%   [matrix column outer, matrix row inner] columns; the reshape/permute
%   below converts between the two.

% DJ (CC), 10/15/2026: Initial coding;

if ~isa(P,'dpvar')
    Ppoly = P;
    return
end
nd = numel(P.dvarname);
m = P.matdim(1);    n = P.matdim(2);
nZ = size(P.degmat,1);
if nd==0
    Ceval = full(P.C);
else
    d0 = zeros(nd,1);
    for k=1:nd
        d0(k) = d0map(P.dvarname{k});
    end
    % Equivalent to left-multiplying by kron(I_m,[1,d0'])
    sel = kron(speye(m),[1,d0(:)']);
    Ceval = full(sel*P.C);
end
% dpvar layout (col = (j-1)*nZ+l) -> 'polynomial' layout (nZ x m*n, col = (j-1)*m+i)
if nZ==0 || m*n==0
    Cpoly = zeros(nZ,m*n);
else
    Cpoly = reshape(permute(reshape(Ceval,[m,nZ,n]),[2,1,3]),[nZ,m*n]);
end
Ppoly = polynomial(Cpoly,P.degmat,P.varname,P.matdim);

end
