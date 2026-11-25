function [nU, ZU, PA, PB, mapA, mapB, unionIsA, unionIsB, dU] = unionBasis(nA, ZA, nB, ZB)
%UNIONBASIS Merge two sorted-unique tensor-product monomial bases.
%
% Convention enforced here:
%   Z_old = P * Z_union
% where P is a sparse row-pick (selection) matrix:
%   - size(P) = d_old x d_union
%   - each row has exactly one 1 selecting the matching union basis element
%
% INPUTS
%   nA, nB : sorted-unique cellstr variable names
%   ZA, ZB : exponent cells aligned to nA, nB, each exponent vector
%            sorted-unique
%
% OUTPUTS
%   nU, ZU        : merged names/exponents (a variable missing from one
%                   operand has exponent [0] there)
%   PA, PB        : sparse selection matrices satisfying:
%                     Z_A = PA * Z_U,   Z_B = PB * Z_U
%   mapA, mapB    : old tensor index -> union tensor index maps (length
%                   dA / dB)
%   unionIsA,     : logical, true iff the union basis equals the basis of
%   unionIsB        A (resp. B), i.e. the map is the identity
%   dU            : union tensor basis size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - unionBasis
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
% SS, 01/16/2026: Initial coding;

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


end

% ---------------- local helpers ----------------

function u = mergeSortedCellstr(a, b)
% Sorted union of two sorted-unique cellstr arrays.
%
% INPUTS
% - a, b:   1 x na and 1 x nb 'cellstr', each sorted and unique;
%
% OUTPUTS
% - u:      1 x k 'cellstr', the sorted union (k <= na+nb).
i = 1; j = 1;
na = numel(a); nb = numel(b);
u = cell(1, na+nb);
k = 0;

while i <= na && j <= nb
    ai = a{i}; bj = b{j};
    % lexicographic order of ai and bj without string allocations: the
    % first differing character decides,
    % else the shorter string is less
    c1 = char(ai); c2 = char(bj); L1 = numel(c1); L2 = numel(c2); L = min(L1, L2);
    dd = c1(1:L) - c2(1:L);
    ix = find(dd ~= 0, 1, 'first');
    if isempty(ix), less = (L1 < L2); else, less = (dd(ix) < 0); end
    if strcmp(ai, bj)
        k=k+1; u{k} = ai; i=i+1; j=j+1;
    elseif less
        k=k+1; u{k} = ai; i=i+1;
    else
        k=k+1; u{k} = bj; j=j+1;
    end
end
while i <= na, k=k+1; u{k} = a{i}; i=i+1; end
while j <= nb, k=k+1; u{k} = b{j}; j=j+1; end

u = u(1:k);
end

function [u, mapA, mapB, unionIsA, unionIsB] = mergeSortedUnique(a, b)
% Merge two sorted-unique numeric vectors a,b into sorted-unique u.
% Returns maps old->union: u(mapA(i)) = a(i), u(mapB(j)) = b(j)
%
% INPUTS
% - a, b:   numeric vectors of na and nb elements, each sorted and unique;
%
% OUTPUTS
% - u:          k x 1 sorted union of a and b;
% - mapA, mapB: na x 1 and nb x 1 positions of the elements of a, b in u;
% - unionIsA, unionIsB: logical, true if u equals a (resp. b).
a = a(:); b = b(:);
na = numel(a); nb = numel(b);

u = zeros(na+nb,1);
mapA = zeros(na,1);
mapB = zeros(nb,1);

% Take the smaller head (both if equal), recording its position in u
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
% Append the remaining tail
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

% The union equals an operand if it has the same length, with identity map
unionIsA = (k==na) && all(mapA == (1:na)');
unionIsB = (k==nb) && all(mapB == (1:nb)');
end

function mapFull = tensorIndexMap(map_cells, nOld_sizes, nU_sizes)
% Build old->union linear index map for tensor-product ordering.
% Ordering matches MATLAB kron: last variable varies fastest.
%
% INPUTS
% - map_cells:  1 x nv cell, map_cells{r} the positions of the old
%               exponents of variable r among its union exponents;
% - nOld_sizes: 1 x nv array, the old number of exponents per variable;
% - nU_sizes:   1 x nv array, the union number of exponents per variable;
%
% OUTPUTS
% - mapFull:    prod(nOld_sizes) x 1 array, the union linear index of each
%               old linear index.

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
