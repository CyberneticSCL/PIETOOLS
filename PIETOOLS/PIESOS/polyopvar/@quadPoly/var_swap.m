function P = var_swap(P, var1, var2)
%VAR_SWAP Swap one left variable's exponent with one right variable's
% exponent in a quadPoly stored with tensor-product monomial bases (cell
% Zs/Zt).
%
% If var1 is not in P.ns, it is appended to P.ns as a constant-only
% variable (degree 0), and likewise var2 to P.nt; an appended name is not
% sorted into place.
%
% After swapping, each exponent list keeps only the degrees that appear,
% and C is remapped accordingly.
%
% INPUTS
% - P:      m x n 'quadPoly' object;
% - var1:   'char', the name of a left variable (in P.ns or new);
% - var2:   'char', the name of a right variable (in P.nt or new);
%
% OUTPUTS
% - P:      m x n 'quadPoly' object with the exponents of var1 and var2
%           swapped.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - var_swap
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
% SS, 01/20/2026: Initial coding;

    ns = P.ns; nt = P.nt;
    Zs = P.Zs; Zt = P.Zt;
    C  = sparse(P.C);

    m = P.dim(1);
    n = P.dim(2);

    % --- normalize empties ---
    if isempty(Zs), Zs = {}; end
    if isempty(Zt), Zt = {}; end
    if isempty(ns), ns = cell(1,numel(Zs)); end
    if isempty(nt), nt = cell(1,numel(Zt)); end

    % basis sizes
    ds = basisSize(Zs);
    dt = basisSize(Zt);

    % ---- ensure var1 exists in ns / Zs, var2 exists in nt / Zt ----
    idxS = find(strcmp(ns, var1), 1);
    if isempty(idxS)
        ns{end+1} = var1;
        Zs{end+1} = 0;           % constant-only degrees
        idxS = numel(ns);
        % update ds
        ds = basisSize(Zs);
    end

    idxT = find(strcmp(nt, var2), 1);
    if isempty(idxT)
        nt{end+1} = var2;
        Zt{end+1} = 0;           % constant-only degrees
        idxT = numel(nt);
        % update dt
        dt = basisSize(Zt);
    end

    % refresh locals after possible append
    P.ns = ns; P.nt = nt;
    P.Zs = Zs; P.Zt = Zt;

    % ---- if C is empty/zero, just commit the variable additions and return ----
    if isempty(C) || nnz(C)==0
        P.C = sparse(m*ds, n*dt);
        return;
    end

    % ---- find all nonzeros ----
    [r_idx, c_idx, val] = find(C);

    % Decode matrix row/col and full monomial indices
    % Decode matrix row/col and full monomial indices
    rIn  = ceil(r_idx / ds);     % 1..m
    cIn  = ceil(c_idx / dt);      % 1..n
    kOld = mod(r_idx - 1, ds) + 1;           % 1..ds (old full s-monomial)
    lOld = mod(c_idx - 1, dt) + 1;           % 1..dt (old full t-monomial)

    % Unique (k,l) pairs
    pairKey = kOld + ds*(lOld-1);
    [uKey, ~, loc] = unique(pairKey);
    ku = mod(uKey - 1, ds) + 1;
    lu = floor((uKey - 1) / ds) + 1;

    % ---- convert full indices to per-variable position tuples ----
    sSizes = cellfun(@numel, Zs);
    tSizes = cellfun(@numel, Zt);

    sSubs = tensorDecomp(ku, sSizes);   % (#pairs)×|ns|
    tSubs = tensorDecomp(lu, tSizes);   % (#pairs)×|nt|

    % ---- convert positions -> actual degrees for each variable ----
    sDeg = subsToDegrees(sSubs, Zs);    % (#pairs)×|ns|
    tDeg = subsToDegrees(tSubs, Zt);    % (#pairs)×|nt|

    % ---- swap degrees between var1(left) and var2(right) ----
    tmp = sDeg(:, idxS);
    sDeg(:, idxS) = tDeg(:, idxT);
    tDeg(:, idxT) = tmp;

    % ---- rebuild per-variable degree lists to include only needed degrees ----
    [Zs_new, sPosMap] = rebuildBasisFromDegrees(Zs, sDeg);
    [Zt_new, tPosMap] = rebuildBasisFromDegrees(Zt, tDeg);

    ds_new = basisSize(Zs_new);
    dt_new = basisSize(Zt_new);

    % ---- map each unique pair to new full indices ----
    % old pair u -> new per-var positions of its degrees in Zs_new/Zt_new
    sSubs_new = degreesToSubs(sDeg, Zs_new);
    tSubs_new = degreesToSubs(tDeg, Zt_new);

    kNew_u = subv2lin_rows(sSubs_new, cellfun(@numel, Zs_new)); % (#pairs)×1 in 1..ds_new
    lNew_u = subv2lin_rows(tSubs_new, cellfun(@numel, Zt_new)); % (#pairs)×1 in 1..dt_new

    % Map each nonzero entry via its pair id
    kNew_nz = kNew_u(loc);
    lNew_nz = lNew_u(loc);

    % ---- rebuild C (sparse sums duplicates automatically) ----

    r_new = kNew_nz + (rIn - 1)*ds_new;
    c_new = lNew_nz + (cIn - 1)*dt_new;

    C_new = sparse(r_new, c_new, val, m*ds_new, n*dt_new);

    % ---- commit ----
    P.ns = ns;
    P.nt = nt;
    P.Zs = Zs_new;
    P.Zt = Zt_new;
    P.C  = C_new;
end

% ====================== helpers ======================

function d = basisSize(Zcell)
% D = BASISSIZE(ZCELL) returns the number of monomials of the tensor
% product basis defined by the exponent vectors in ZCELL (an empty
% element counts as the constant monomial).
%
% INPUTS
% - Zcell:  1 x k cell of exponent vectors;
%
% OUTPUTS
% - d:      positive integer, the product of their lengths (1 if empty).
    if isempty(Zcell), d = 1; return; end
    d = 1;
    for i = 1:numel(Zcell)
        zi = Zcell{i};
        if isempty(zi), zi = 0; end
        d = d * numel(zi);
    end
end

function subs = tensorDecomp(idx, sizes)
% SUBS = TENSORDECOMP(IDX,SIZES) returns the subscripts of the linear
% indices IDX in an array of size SIZES, the last dimension varying
% fastest.
%
% INPUTS
% - idx:    vector of linear indices in 1..prod(sizes);
% - sizes:  1 x k array of positive integers;
%
% OUTPUTS
% - subs:   length(idx) x k array of subscripts (1-based).
    idx = idx(:) - 1;
    k = numel(sizes);
    subs = zeros(numel(idx), k);
    for i = k:-1:1
        si = sizes(i);
        subs(:,i) = mod(idx, si) + 1;
        idx = floor(idx / si);
    end
end

function deg = subsToDegrees(subs, Zcell)
% DEG = SUBSTODEGREES(SUBS,ZCELL) returns the degrees at the positions SUBS
% in the exponent lists ZCELL (an empty list counts as the degree 0).
%
% INPUTS
% - subs:   L x k array of positions;
% - Zcell:  1 x k cell of exponent vectors;
%
% OUTPUTS
% - deg:    L x k array of degrees.
    L = size(subs,1);
    k = size(subs,2);
    deg = zeros(L,k);
    for i = 1:k
        zi = Zcell{i}(:);
        if isempty(zi), zi = 0; end
        deg(:,i) = zi(subs(:,i));
    end
end

function [Znew, posMap] = rebuildBasisFromDegrees(Zold, degUsed)
% [ZNEW,POSMAP] = REBUILDBASISFROMDEGREES(ZOLD,DEGUSED) sets, for each
% variable i, Znew{i} = unique(degUsed(:,i),'sorted') (0 if empty).
%
% INPUTS
% - Zold:     1 x k cell of exponent vectors (only its length is used);
% - degUsed:  L x k array of the degrees used;
%
% OUTPUTS
% - Znew:     1 x k cell of sorted unique degree vectors;
% - posMap:   1 x k cell of empty arrays (not filled).
    k = numel(Zold);
    Znew = cell(1,k);
    posMap = cell(1,k); %#ok<NASGU>
    for i = 1:k
        z = unique(degUsed(:,i), 'sorted');
        if isempty(z), z = 0; end
        Znew{i} = z;
    end
end

function subs = degreesToSubs(deg, Zcell)
% SUBS = DEGREESTOSUBS(DEG,ZCELL) returns the positions of the degrees DEG
% in the exponent lists ZCELL; a degree missing from its list is an error.
%
% INPUTS
% - deg:    L x k array of degrees;
% - Zcell:  1 x k cell of exponent vectors;
%
% OUTPUTS
% - subs:   L x k array of positions in the lists of Zcell.
    L = size(deg,1);
    k = size(deg,2);
    subs = zeros(L,k);
    for i = 1:k
        zi = Zcell{i}(:);
        % degrees were built from zi so ismember should succeed
        [tf, loc] = ismember(deg(:,i), zi);
        if any(~tf)
            error('var_swap:degreeMissing', 'A needed degree is missing in rebuilt basis.');
        end
        subs(:,i) = loc;
    end
end

function idx = subv2lin_rows(subs, sizes)
% IDX = SUBV2LIN_ROWS(SUBS,SIZES) returns the linear index of each row of
% subscripts SUBS in an array of size SIZES, the last dimension varying
% fastest.
%
% INPUTS
% - subs:   L x k array of subscripts (1-based);
% - sizes:  1 x k array of positive integers;
%
% OUTPUTS
% - idx:    L x 1 array of linear indices in 1..prod(sizes).
    L = size(subs,1);
    k = numel(sizes);
    idx0 = zeros(L,1);
    stride = ones(1,k);
    for i = 1:k-1
        stride(i) = prod(sizes(i+1:end));
    end
    stride(k) = 1;

    for i = 1:k
        idx0 = idx0 + (subs(:,i)-1) * stride(i);
    end
    idx = idx0 + 1;
end
