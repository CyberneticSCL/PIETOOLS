function H = horzcat(varargin)
%horzcat Horizontal concatenation [A B C ...] for quadPoly.
%
% Inputs:
%   varargin: cell structure of polynomials to be concatenated {A,B,C...}
%
% Outputs:
%   H: quadPoly object representing [A,B,C...]
% Rules:
%   - All inputs must have the same number of rows (matrix dimension m).
%   - Left variables (ns/Zs) stay on the left; right variables (nt/Zt) stay on the right.
%   - Bases are merged (union) across all inputs before concatenation.
%   - Uses unionBasis convention: Z_old = P * Z_union (row-pick), but we lift via maps.
%
% Supports mixing numeric matrices with quadPoly (numeric treated as constant quadPoly).


% ---- Convert all inputs to quadPoly ----
Q = cell(size(varargin));
for k = 1:nargin
    if isa(varargin{k},'quadPoly')
        Q{k} = varargin{k}
    elseif isnumeric(X)
        A = sparse(varargin{k});
        Q{k} = quadPoly(A, {}, {}, size(A), {}, {});
    else
        error('quadPoly:horzcat:badType', 'horzcat supports quadPoly and numeric inputs.');
    end
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
        Ccat = liftIndex(Ccat, m, nCat, mapS_cat, mapT_cat, dsNew, dtNew);
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
        CP = liftIndex(P.C, m, nP, mapS_P, mapT_P, dsU, dtU);
    end

    % Concatenate along columns (same m, same dsU/dtU after lifting)
    Ccat = [Ccat, CP]; %#ok<AGROW>
    nCat = nCat + nP;
end

H = quadPoly(Ccat, ZsU, ZtU, [m, nCat], nsU, ntU);
end
