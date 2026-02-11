function [nU, ZU, PA, PB, mapA, mapB, unionIsA, unionIsB, dU] = unionBasis(nA, ZA, nB, ZB)
%UNIONBASIS Merge two tensor-product monomial bases (sorted-unique names/exponents).
%
% Convention enforced here:
%   Z_old = P * Z_union
% where P is a sparse row-pick (selection) matrix:
%   - size(P) = d_old × d_union
%   - each row has exactly one 1 selecting the matching union basis element
%
% Inputs:
%   nA, nB : sorted-unique cellstr variable names
%   ZA, ZB : exponent cells aligned to nA,nB, each exponent vector sorted-unique
%
% Outputs:
%   nU, ZU        : merged names/exponents (missing variable => exponent [0])
%   PA, PB        : sparse selection matrices satisfying:
%                     Z_A = PA * Z_U,   Z_B = PB * Z_U
%   mapA, mapB    : old tensor index -> union tensor index maps (length dA / dB)
%   unionIsA/B    : true iff union basis equals operand basis (identity map)
%   dU            : union tensor basis size

% ---- merge variable names (sorted union) ----
nU = mergeSortedCellstr(nA, nB);
nv = numel(nU);

ZU = cell(1, nv);

% per-variable maps and sizes (to build full tensor maps without kron)
mapA_cells = cell(1, nv);
mapB_cells = cell(1, nv);
nA_sizes   = zeros(1, nv);
nB_sizes   = zeros(1, nv);
nU_sizes   = zeros(1, nv);

unionIsA = true;
unionIsB = true;

% pointers since nA,nB are sorted
iA = 1; iB = 1;
nNA = numel(nA);
nNB = numel(nB);

for k = 1:nv
    v = nU{k};

    % exponents for v in A (or [0] if missing)
    if iA <= nNA && all(strcmp(nA{iA}, v))
        zA = ZA{iA}(:); iA = iA + 1;
    else
        zA = 0;
    end

    % exponents for v in B (or [0] if missing)
    if iB <= nNB && all(strcmp(nB{iB}, v))
        zB = ZB{iB}(:); iB = iB + 1;
    else
        zB = 0;
    end

    % merge per-variable exponents + maps old->union
    [zU, mA, mB, isA, isB] = mergeSortedUnique(zA, zB);

    ZU{k} = zU;
    mapA_cells{k} = mA;  % length numel(zA), values 1..numel(zU)
    mapB_cells{k} = mB;

    nA_sizes(k) = numel(zA);
    nB_sizes(k) = numel(zB);
    nU_sizes(k) = numel(zU);

    unionIsA = unionIsA && isA;
    unionIsB = unionIsB && isB;
end

dU = prod(nU_sizes);

% Build full tensor old->union maps (length dA / dB)
if unionIsA
    mapA = (1:dU).';
else
    mapA = tensorIndexMap(mapA_cells, nA_sizes, nU_sizes);
end

if unionIsB
    mapB = (1:dU).';
else
    mapB = tensorIndexMap(mapB_cells, nB_sizes, nU_sizes);
end

% Build sparse selection matrices satisfying Z_old = P * Z_union
% (row-pick: row r picks column map(r))
dA = prod(nA_sizes);
dB = prod(nB_sizes);
PA = sparse(1:dA, mapA, 1, dA, dU);
PB = sparse(1:dB, mapB, 1, dB, dU);


% maybe we do not need PA and PB, ever?
end

% ---------------- local helpers ----------------

function u = mergeSortedCellstr(a, b)
% Sorted union of two sorted-unique cellstr arrays.
i = 1; j = 1;
na = numel(a); nb = numel(b);
u = cell(1, na+nb);
k = 0;

while i <= na && j <= nb
    ai = a{i}; bj = b{j};
    if strcmp(ai, bj)
        k=k+1; u{k} = ai; i=i+1; j=j+1;
    elseif lexLessChar(ai, bj)
        k=k+1; u{k} = ai; i=i+1;
    else
        k=k+1; u{k} = bj; j=j+1;
    end
end
while i <= na, k=k+1; u{k} = a{i}; i=i+1; end
while j <= nb, k=k+1; u{k} = b{j}; j=j+1; end

u = u(1:k);
end

function tf = lexLessChar(s1, s2)
% Lexicographic compare without string allocations.
s1 = char(s1); s2 = char(s2);
L1 = numel(s1); L2 = numel(s2);
L  = min(L1, L2);

d = s1(1:L) - s2(1:L);
idx = find(d ~= 0, 1, 'first');
if isempty(idx)
    tf = (L1 < L2);
else
    tf = (d(idx) < 0);
end
end

function [u, mapA, mapB, unionIsA, unionIsB] = mergeSortedUnique(a, b)
% Merge two sorted-unique numeric vectors a,b into sorted-unique u.
% Returns maps old->union: u(mapA(i)) = a(i), u(mapB(j)) = b(j)
a = a(:); b = b(:);
na = numel(a); nb = numel(b);

u = zeros(na+nb,1);
mapA = zeros(na,1);
mapB = zeros(nb,1);

i=1; j=1; k=0;
while i<=na && j<=nb
    if a(i) == b(j)
        k=k+1; u(k)=a(i);
        mapA(i)=k; mapB(j)=k;
        i=i+1; j=j+1;
    elseif a(i) < b(j)
        k=k+1; u(k)=a(i);
        mapA(i)=k;
        i=i+1;
    else
        k=k+1; u(k)=b(j);
        mapB(j)=k;
        j=j+1;
    end
end
while i<=na
    k=k+1; u(k)=a(i);
    mapA(i)=k;
    i=i+1;
end
while j<=nb
    k=k+1; u(k)=b(j);
    mapB(j)=k;
    j=j+1;
end

u = u(1:k);

unionIsA = (k==na) && all(mapA == (1:na)');
unionIsB = (k==nb) && all(mapB == (1:nb)');
end

function mapFull = tensorIndexMap(map_cells, nOld_sizes, nU_sizes)
% Build old->union linear index map for tensor-product ordering.
% Ordering matches MATLAB kron: last variable varies fastest.

nv = numel(map_cells);

mapFull     = 1;   % mapping for empty tail
tailU       = 1;   % product of union sizes in tail
tailOld     = 1;   % product of old sizes in tail

for r = nv:-1:1
    map_r  = map_cells{r}(:);        % length nOld_r, values 1..nU_r
    nOld_r = nOld_sizes(r);

    % idxU = (map_r-1)*tailU + idxTail
    mapFull = kron(ones(nOld_r,1), mapFull) + tailU * kron(map_r-1, ones(tailOld,1));

    tailU   = nU_sizes(r) * tailU;
    tailOld = nOld_r * tailOld;
end

mapFull = mapFull(:);
end
