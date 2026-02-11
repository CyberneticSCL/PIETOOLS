function H = vertcat(varargin)
%vertcat Vertical concatenation [A; B; C; ...] for quadPoly.
%
% Inputs:
%   varargin: cell structure of polynomials to be concatenated {A,B,C...}
%
% Outputs:
%   H: quadPoly object representing [A;B;C...]
% Rules:
%   - All inputs must have the same number of columns (matrix dimension n).
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
        error('quadPoly:vertcat:badType', 'vertcat supports quadPoly and numeric inputs.');
    end
end

% ---- Dimension check (same #cols) ----
n = Q{1}.dim(2);
for k = 2:numel(Q)
    if Q{k}.dim(2) ~= n
        error('quadPoly:vertcat:dimMismatch', 'All inputs must have same number of columns.');
    end
end

% ---- Initialize accumulator with first block ----
nsU = Q{1}.ns;  ZsU = Q{1}.Zs;
ntU = Q{1}.nt;  ZtU = Q{1}.Zt;

Ccat = sparse(Q{1}.C);
mCat = Q{1}.dim(1);

% ---- Fold-in remaining blocks, updating union bases and lifting as needed ----
for k = 2:numel(Q)
    P = Q{k};
    mP = P.dim(1);

    % --- merge s-side bases: current union with P ---
    [nsNew, ZsNew, ~, ~, mapS_cat, mapS_P, sIsCat, sIsP, dsNew] = unionBasis(nsU, ZsU, P.ns, P.Zs);

    % --- merge t-side bases: current union with P ---
    [ntNew, ZtNew, ~, ~, mapT_cat, mapT_P, tIsCat, tIsP, dtNew] = unionBasis(ntU, ZtU, P.nt, P.Zt);

    % If union changed, lift current concatenated coefficients to new union basis
    if ~(sIsCat && tIsCat)
        Ccat = liftIndex(Ccat, mCat, n, mapS_cat, mapT_cat, dsNew, dtNew);
        nsU = nsNew; ZsU = ZsNew; dsU = dsNew;
        ntU = ntNew; ZtU = ZtNew; dtU = dtNew;
    else
        nsU = nsNew; ZsU = ZsNew; dsU = dsNew;
        ntU = ntNew; ZtU = ZtNew; dtU = dtNew;
    end

    % Lift P into current union basis (if needed)
    if sIsP && tIsP
        CP = sparse(P.C);
    else
        CP = liftIndex(P.C, mP, n, mapS_P, mapT_P, dsU, dtU);
    end

    % Concatenate along rows (same n, same dsU/dtU after lifting)
    Ccat = [Ccat; CP]; %#ok<AGROW>
    mCat = mCat + mP;
end

H = quadPoly(Ccat, ZsU, ZtU, [mCat, n], nsU, ntU);
end

