function Pnew = lmult_params(M,P,nd)
% PNEW = LMULT_PARAMS(M,P,ND) computes the coefficient matrix of the
% product of the mm x p matrix M with the p x n matrix of polynomials
% defined by the sparse coefficient matrix P, with rows (i-1)*(nd+1)+d for
% matrix row i and decision variable d (constant term first), as for
% 'dpvar'.
%
% INPUTS
% - M:      mm x p array of type 'double' (sparse or full);
% - P:      p*(nd+1) x c sparse array defining the p x n matrix of
%           polynomials. The columns of P are not mixed by the product;
% - nd:     scalar integer specifying the number of decision variables;
%
% OUTPUTS
% - Pnew:   mm*(nd+1) x c sparse array defining the mm x n matrix of
%           polynomials M*P, with the rows ordered as for P;
%

% DJ, 10/15/2026: Initial coding;

[mm,p] = size(M);
M = sparse(M);
[i,j,v] = find(P);
i = i(:);   j = j(:);   v = v(:);       % row vectors if P has one row
blk = ceil(i/(nd+1));
d = i - (blk-1)*(nd+1);
if p==1
    ord = (1:numel(i))';        edges = [0;numel(i)];
else
    [~,ord] = sort(blk);
    edges = [0;cumsum(accumarray(blk,1,[p,1]))];
end
[ia,ja,va] = find(M);
Ii = cell(numel(ia),1);     Jj = Ii;    Vv = Ii;
for t=1:numel(ia)
    idx = ord(edges(ja(t))+1:edges(ja(t)+1));
    Ii{t} = (ia(t)-1)*(nd+1) + d(idx);
    Jj{t} = j(idx);
    Vv{t} = v(idx)*va(t);
end
Pnew = sparse(vertcat(Ii{:}),vertcat(Jj{:}),vertcat(Vv{:}),mm*(nd+1),size(P,2));

end
