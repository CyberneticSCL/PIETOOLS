function H = horzcat(varargin)
%HORZCAT Horizontal concatenation [A B C ...] for quadPoly.
%
% Rules:
%   - All inputs must have the same number of rows (matrix dimension m).
%   - Left variables (ns/Zs) stay on the left; right variables (nt/Zt) stay on the right.
%   - Bases are merged (union) across all inputs before concatenation.
%   - Uses unionBasis convention: Z_old = P * Z_union (row-pick), but we lift via maps.
%
% Supports mixing numeric matrices with quadPoly (numeric treated as constant quadPoly).

if nargin == 0
    H = quadPoly(sparse([]), {}, {}, [0 0], {}, {});
    return;
end

% ---- Convert all inputs to quadPoly ----
Q = cell(size(varargin));
for k = 1:nargin
    Q{k} = toQuadPoly(varargin{k});
end

% ---- Dimension check (same #rows) ----
m = Q{1}.dim(1);
for k = 2:numel(Q)
    if Q{k}.dim(1) ~= m
        error('quadPoly:horzcat:dimMismatch', 'All inputs must have same number of rows.');
    end
end

% ---- Initialize accumulator with first block ----
nsU = Q{1}.ns;  ZsU = Q{1}.Zs;
ntU = Q{1}.nt;  ZtU = Q{1}.Zt;

dsU = prod(cellfun(@numel, ZsU));
dtU = prod(cellfun(@numel, ZtU));

Ccat = sparse(Q{1}.C);
nCat = Q{1}.dim(2);

% ---- Fold-in remaining blocks, updating union bases and lifting as needed ----
for k = 2:numel(Q)
    P = Q{k};
    nP = P.dim(2);

    % --- merge s-side bases: current union with P ---
    [nsNew, ZsNew, ~, ~, mapS_cat, mapS_P, sIsCat, sIsP, dsNew] = unionBasis(nsU, ZsU, P.ns, P.Zs);

    % --- merge t-side bases: current union with P ---
    [ntNew, ZtNew, ~, ~, mapT_cat, mapT_P, tIsCat, tIsP, dtNew] = unionBasis(ntU, ZtU, P.nt, P.Zt);

    % If union changed, lift current concatenated coefficients to new union basis
    if ~(sIsCat && tIsCat)
        Ccat = liftIndexMaps(Ccat, m, nCat, mapS_cat, mapT_cat, dsNew, dtNew);
        nsU = nsNew; ZsU = ZsNew; dsU = dsNew;
        ntU = ntNew; ZtU = ZtNew; dtU = dtNew;
    else
        % union unchanged (still update names/exps for consistency)
        nsU = nsNew; ZsU = ZsNew; dsU = dsNew;
        ntU = ntNew; ZtU = ZtNew; dtU = dtNew;
    end

    % Lift P into current union basis (if needed)
    if sIsP && tIsP
        CP = sparse(P.C);
    else
        CP = liftIndexMaps(P.C, m, nP, mapS_P, mapT_P, dsU, dtU);
    end

    % Concatenate along columns (same m, same dsU/dtU after lifting)
    Ccat = [Ccat, CP]; %#ok<AGROW>
    nCat = nCat + nP;
end

H = quadPoly(Ccat, ZsU, ZtU, [m, nCat], nsU, ntU);

end

% =========================================================================
% Local helpers
% =========================================================================
function Q = toQuadPoly(X)
% Convert numeric to constant quadPoly; pass-through quadPoly.
if isa(X,'quadPoly')
    Q = X;
elseif isnumeric(X)
    A = sparse(X);
    Q = quadPoly(A, {}, {}, size(A), {}, {});
else
    error('quadPoly:horzcat:badType', 'horzcat supports quadPoly and numeric inputs.');
end
end

function Cup = liftIndexMaps(C, m, n, mapS, mapT, dsU, dtU)
% Lift coefficient matrix by reindexing nnz entries using old->union maps.
C = sparse(C);

dsOld = numel(mapS);
dtOld = numel(mapT);

[i, j, v] = find(C);

a  = mod(i-1, dsOld) + 1;
ib = floor((i-1) / dsOld);

b  = mod(j-1, dtOld) + 1;
jb = floor((j-1) / dtOld);

i2 = ib * dsU + mapS(a);
j2 = jb * dtU + mapT(b);

Cup = sparse(i2, j2, v, m*dsU, n*dtU);
end
