function H = mtimes(F, G)
%MTIMES Matrix product of quadPoly objects (and/or numeric matrices).
%
% Convention for basis selection matrices throughout the codebase:
%   Z_old = P * Z_union
%
% For quadPoly * quadPoly:
%   - we first align (union + lift) the s- and t- bases so both operands share them
%   - then we multiply coefficients with exponent-addition (Minkowski sum) on each
%     tensor-decomposed variable, i.e., exponents add under multiplication.
%   - IMPORTANT: left variables stay left (ns/Zs), right variables stay right (nt/Zt).
%
% Numeric cases:
%   A*quadPoly or quadPoly*B are supported (A,B numeric matrices).

% ---------- numeric dispatch ----------
if isnumeric(F) && isa(G,'quadPoly')
    H = leftNumericTimes(F, G);
    return;
elseif isa(F,'quadPoly') && isnumeric(G)
    H = rightNumericTimes(F, G);
    return;
elseif ~(isa(F,'quadPoly') && isa(G,'quadPoly'))
    error('quadPoly:mtimes:badType', 'mtimes requires quadPoly or numeric.');
end

% ---------- quadPoly * quadPoly ----------
m = F.dim(1);
n = F.dim(2);
n2 = G.dim(1);
p = G.dim(2);
if n ~= n2
    error('quadPoly:mtimes:dimMismatch', 'Inner dimensions must match.');
end

% --- 1) Align variable sets + bases (s-side and t-side) ---
[nsU, ZsU, ~, ~, mapSF, mapSG, sIsF, sIsG, dsU] = unionBasis(F.ns, F.Zs, G.ns, G.Zs);
[ntU, ZtU, ~, ~, mapTF, mapTG, tIsF, tIsG, dtU] = unionBasis(F.nt, F.Zt, G.nt, G.Zt);

if ~(isequal(nsU, nsU) && isequal(ntU, ntU)) %#ok<*ISEQ>
    % kept for clarity; nsU/ntU are the aligned variable lists
end

% Lift coefficients into the aligned (union) bases
if sIsF && tIsF
    CF = F.C;
else
    CF = liftIndexMaps(F.C, m, n, mapSF, mapTF, dsU, dtU);
end
if sIsG && tIsG
    CG = G.C;
else
    CG = liftIndexMaps(G.C, n, p, mapSG, mapTG, dsU, dtU);
end

% After lifting:
%   CF is (m*dsU) × (n*dtU)
%   CG is (n*dsU) × (p*dtU)

% --- 2) Build product bases (exponents add under multiplication) ---
% For each variable, product exponents are the sorted-unique Minkowski sum.
[ZsP, sMapMats, dsP, sStrideP] = productBasis1D(ZsU, ZsU);  % left variables: add exponents
[ZtP, tMapMats, dtP, tStrideP] = productBasis1D(ZtU, ZtU);  % right variables

% --- 3) Multiply blocks with exponent-addition mapping ---
% H will be (m*dsP) × (p*dtP)
CH = spalloc(m*dsP, p*dtP, nnz(CF) + nnz(CG));  % rough initial guess

% Precompute sizes per variable for index decomposition
sSizes = cellfun(@numel, ZsU);  % same for both operands after union
tSizes = cellfun(@numel, ZtU);

for i = 1:m
    rF = (i-1)*dsU + (1:dsU);               % rows for block-row i in CF

    for j = 1:p
        % We will accumulate the dsP×dtP block for (i,j)
        blockRows = (i-1)*dsP;
        blockCols = (j-1)*dtP;

        % Accumulate contributions over k
        Iacc = [];
        Jacc = [];
        Vacc = [];

        for k = 1:n
            cF = (k-1)*dtU + (1:dtU);
            rG = (k-1)*dsU + (1:dsU);
            cG = (j-1)*dtU + (1:dtU);

            BF = CF(rF, cF);  % dsU×dtU
            BG = CG(rG, cG);  % dsU×dtU

            if nnz(BF) == 0 || nnz(BG) == 0
                continue;
            end

            [a, b, v1] = find(BF);   % indices in 1..dsU and 1..dtU
            [c, d, v2] = find(BG);

            % Form all pairs of terms BF(a,b)*BG(c,d)
            % Pair arrays length = nnz(BF)*nnz(BG)
            na = numel(v1); nc = numel(v2);

            aP = repelem(a, nc);
            bP = repelem(b, nc);
            vP = repelem(v1, nc) .* repmat(v2, na, 1);

            cP = repmat(c, na, 1);
            dP = repmat(d, na, 1);

            % Map (a,c) -> s-basis product index, (b,d) -> t-basis product index
            sIdx = pairToProdIndex(aP, cP, sSizes, sMapMats, sStrideP, dsP);
            tIdx = pairToProdIndex(bP, dP, tSizes, tMapMats, tStrideP, dtP);

            Iacc = [Iacc; blockRows + sIdx]; %#ok<AGROW>
            Jacc = [Jacc; blockCols + tIdx]; %#ok<AGROW>
            Vacc = [Vacc; vP];               %#ok<AGROW>
        end

        if ~isempty(Vacc)
            CH = CH + sparse(Iacc, Jacc, Vacc, m*dsP, p*dtP);
        end
    end
end

H = quadPoly(CH, ZsP, ZtP, [m p], nsU, ntU);

end

% =========================================================================
% Numeric helpers
% =========================================================================
function H = leftNumericTimes(A, F)
% A is (q×m), F is (m×n) quadPoly -> H is (q×n) quadPoly
[q, m] = size(A);
if m ~= F.dim(1)
    error('quadPoly:mtimes:leftNumericDimMismatch', 'Left numeric matrix has wrong size.');
end

ds = prod(cellfun(@numel, F.Zs));
dt = prod(cellfun(@numel, F.Zt));

% Left multiplication: (A*F)(s,t) = A * F(s,t)
% Coeff update: Cnew = (kron(A, I_ds)) * C
L = kron(sparse(A), speye(ds));
Cnew = L * sparse(F.C);

H = quadPoly(Cnew, F.Zs, F.Zt, [q, F.dim(2)], F.ns, F.nt);
end

function H = rightNumericTimes(F, B)
% F is (m×n) quadPoly, B is (n×p) numeric -> H is (m×p) quadPoly
[n, p] = size(B);
if n ~= F.dim(2)
    error('quadPoly:mtimes:rightNumericDimMismatch', 'Right numeric matrix has wrong size.');
end

ds = prod(cellfun(@numel, F.Zs));
dt = prod(cellfun(@numel, F.Zt));

% Right multiplication: (F*B)(s,t) = F(s,t) * B
% Coeff update: Cnew = C * kron(B, I_dt)
R = kron(sparse(B), speye(dt));
Cnew = sparse(F.C) * R;

H = quadPoly(Cnew, F.Zs, F.Zt, [F.dim(1), p], F.ns, F.nt);
end

% =========================================================================
% Basis product helpers (tensor-decomposed)
% =========================================================================
function [ZP, mapMats, dP, strideP] = productBasis1D(ZA, ZB)
%PRODUCTBASIS1D Per-variable Minkowski sums and lookup maps for tensor basis product.
%
% Inputs:
%   ZA, ZB : 1×k cells, sorted-unique exponent vectors per variable
% Outputs:
%   ZP      : 1×k cells, sorted-unique sums per variable
%   mapMats : 1×k cells, mapMats{i} is (nAi×nBi) matrix giving index in ZP{i}
%   dP      : total tensor basis size prod_i numel(ZP{i})
%   strideP : 1×k strides for linear indexing in kron order (last varies fastest)

k = numel(ZA);
ZP = cell(1,k);
mapMats = cell(1,k);

sizesP = zeros(1,k);

for i = 1:k
    a = ZA{i}(:);
    b = ZB{i}(:);

    % All pairwise sums in this variable
    S = a + b.';                         % nAi × nBi
    z = unique(S(:), 'sorted');          % sorted-unique sums
    ZP{i} = z;
    sizesP(i) = numel(z);

    % Map each sum back to its index in z
    [~, pos] = ismember(S(:), z);
    mapMats{i} = reshape(pos, size(S,1), size(S,2));
end

dP = prod(sizesP);

% Strides for kron order (last variable varies fastest)
strideP = ones(1,k);
for i = 1:k-1
    strideP(i) = prod(sizesP(i+1:end));
end
strideP(k) = 1;
end

function idxP = pairToProdIndex(idxA, idxB, sizes, mapMats, strideP, dP)
%PAIRTOPRODINDEX Map pairs of tensor indices (idxA,idxB) -> product tensor index.
%
% idxA, idxB are vectors (same length L) with values in 1..prod(sizes).
% sizes   : 1×k sizes of the aligned tensor basis per variable
% mapMats : 1×k, mapMats{i}(a_i,b_i) gives per-variable index in product basis
% strideP : 1×k strides for linear index in product basis (kron order)
%
% Returns idxP in 1..dP.

L = numel(idxA);
k = numel(sizes);

subsA = tensorDecomp(idxA, sizes);  % L×k (each entry 1..sizes(i))
subsB = tensorDecomp(idxB, sizes);

idx0 = zeros(L,1);
for i = 1:k
    map_i = mapMats{i};
    si = map_i(subsA(:,i), subsB(:,i));  % L×1 indices in 1..numel(ZP{i})
    idx0 = idx0 + (si - 1) * strideP(i);
end

idxP = idx0 + 1;
% (Optional) trust invariants; idxP should be within 1..dP
end

function subs = tensorDecomp(idx, sizes)
%TENSORDECOMP Convert linear tensor index -> per-variable subscripts (kron order).
%
% Ordering matches MATLAB kron: last variable varies fastest.
% idx is 1-based vector.

idx = idx(:) - 1;   % 0-based
k = numel(sizes);
subs = zeros(numel(idx), k);

for i = k:-1:1
    si = sizes(i);
    subs(:,i) = mod(idx, si) + 1;
    idx = floor(idx / si);
end
end

% =========================================================================
% Requires: liftIndexMaps (from earlier)
% =========================================================================
function Cup = liftIndexMaps(C, m, n, mapS, mapT, dsU, dtU)
%LIFTINDEXMAPS Lift coefficient matrix by reindexing nnz entries using maps.
%
% mapS: dsOld×1 old->union map (values in 1..dsU)
% mapT: dtOld×1 old->union map (values in 1..dtU)

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
