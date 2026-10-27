function H = subs_OLD(F, vars, vals)
%SUBS Substitute variables in a quadPoly with tensor-decomposed exponent bases.
%
% Supports substitutions on BOTH sides:
%   s-variables: names in F.ns (typically like 's1','s2',...)
%   t-variables: names in F.nt (typically like 't1','theta2','t3',...)
%
% RHS values may be:
%   - numeric/logical 0 or 1
%   - variable name (string/char) from either side, e.g. 's2','t3','theta2'
%
% Semantics:
%   - Setting a variable to 0 kills any term where its exponent is nonzero.
%   - Setting to 1 drops that exponent.
%   - Mapping to another variable transfers exponent (adds into target exponent).
%
% IMPORTANT: This adapts the older matrix-exponent version to the current representation:
%   F.Zs and F.Zt are 1×k cell arrays, one exponent list per variable.
%   A monomial basis index corresponds to a tensor-product choice of one exponent
%   from each variable's exponent list (last variable varies fastest).

% Normalize input pairs
if nargin < 2 || isempty(vars)
    H = F;
    return;
end
[vars, vals] = normalizePairs(vars, vals);

% Build name->position maps for left and right variables
[nsMap, ks] = nameIndexMap(F.ns);
[ntMap, kt] = nameIndexMap(F.nt);

% Build substitution plans (per original variable)
[modeS,toS,modeT,toT] = buildMaps(ks, kt, nsMap, ntMap, vars, vals);

m = F.dim(1);
n = F.dim(2);

ds = prod(cellfun(@numel, F.Zs));
dt = prod(cellfun(@numel, F.Zt));

[I,J,V] = find(F.C);
if isempty(V)
    % Return zero with trivial bases
    H = quadPoly(sparse(m*1, n*1), {0}, {0}, F.dim, F.ns, F.nt);
    return;
end

% Block row/col indices and local tensor indices in s/t bases
br = floor((I-1)/ds) + 1;   ls = I - (br-1)*ds;  % ls in 1..ds
bc = floor((J-1)/dt) + 1;   lt = J - (bc-1)*dt;  % lt in 1..dt

% We will build NEW tensor-decomposed bases (cells) by collecting all exponents
% that appear after substitution.
ZsNewSets = cell(1, ks);   % each is a vector (will store unique values)
ZtNewSets = cell(1, kt);

% For speed: use MATLAB containers.Map with numeric keys per variable
mapS = cell(1, ks);  % mapS{i}: exponent -> new local index in ZsNew{i}
mapT = cell(1, kt);  % mapT{i}: exponent -> new local index in ZtNew{i}
for i = 1:ks
    mapS{i} = containers.Map('KeyType','double','ValueType','double');
end
for i = 1:kt
    mapT{i} = containers.Map('KeyType','double','ValueType','double');
end

% We also need to map old tensor indices (ls,lt) to new tensor indices after substitution.
% We will:
%   1) decode ls -> exponent vector (per s variable)
%   2) decode lt -> exponent vector (per t variable)
%   3) apply substitution -> new exponent vectors (s and t)
%   4) insert each per-variable exponent into ZsNewSets/ZtNewSets and get local indices
%   5) compose new tensor linear index from those local indices
%
% Cache decoded exponent vectors for repeated ls/lt
cacheS = containers.Map('KeyType','double','ValueType','any'); % ls -> exponents (ks×1)
cacheT = containers.Map('KeyType','double','ValueType','any'); % lt -> exponents (kt×1)

% Accumulate lifted triplets for new coefficient matrix
brH = zeros(numel(V),1);
bcH = zeros(numel(V),1);
isH = zeros(numel(V),1);   % new s tensor index (1..dsN)
itH = zeros(numel(V),1);   % new t tensor index (1..dtN)
vH  = zeros(numel(V),1);

kout = 0;

for k = 1:numel(V)
    % Decode old tensor indices to exponent vectors
    es = getCachedExp(cacheS, ls(k), F.Zs);
    et = getCachedExp(cacheT, lt(k), F.Zt);

    % Apply substitution maps -> new exponent vectors
    [es2, et2, ok] = applyMapsCell(es, et, modeS, toS, modeT, toT);
    if ~ok
        continue;
    end

    % Insert per-variable exponents into new sets and get per-variable local indices
    [subS, ZsNewSets, mapS] = internExpVec(es2, ZsNewSets, mapS);
    [subT, ZtNewSets, mapT] = internExpVec(et2, ZtNewSets, mapT);

    kout = kout + 1;
    brH(kout) = br(k);
    bcH(kout) = bc(k);
    isH(kout) = sub2linTensor(subS, sizesFromSets(ZsNewSets));
    itH(kout) = sub2linTensor(subT, sizesFromSets(ZtNewSets));
    vH(kout)  = V(k);
end

if kout == 0
    % Everything killed
    ZsNew = cellfun(@(~) 0, F.Zs, 'UniformOutput', false);
    ZtNew = cellfun(@(~) 0, F.Zt, 'UniformOutput', false);
    dsN = prod(cellfun(@numel, ZsNew));
    dtN = prod(cellfun(@numel, ZtNew));
    CH  = sparse(m*dsN, n*dtN);
    H   = quadPoly(CH, ZsNew, ZtNew, F.dim, F.ns, F.nt);
    return;
end

% Trim
brH = brH(1:kout);
bcH = bcH(1:kout);
isH = isH(1:kout);
itH = itH(1:kout);
vH  = vH(1:kout);

% Finalize bases: sort unique exponents per variable and rebuild index maps
ZsNew = finalizeSets(ZsNewSets);
ZtNew = finalizeSets(ZtNewSets);

dsN = prod(cellfun(@numel, ZsNew));
dtN = prod(cellfun(@numel, ZtNew));

% Recompute tensor indices isH/itH because finalizeSets sorts (indices may change)
% Build per-variable maps from exponent value -> sorted position
val2idxS = valueToIndexMaps(ZsNew);
val2idxT = valueToIndexMaps(ZtNew);

% Second pass: rebuild new tensor indices consistently with sorted bases
for kk = 1:kout
    es2 = getCachedExp(cacheS, ls(findIndex(kk, kout, ls, br, ds)), F.Zs); %#ok<NASGU>
end
% Instead of trying to reuse old loop indices, rebuild using stored subS/subT would be ideal.
% To keep this function simple and correct, we store the exponents during the first pass.

% ---- rebuild properly by storing exponents in first pass (minimal extra memory) ----
% (We implement that now by redoing first pass with storage.)

% Re-run in a clean way: store substituted exponents for each kept term
esStore = zeros(kout, ks);
etStore = zeros(kout, kt);

kout2 = 0;
cacheS2 = containers.Map('KeyType','double','ValueType','any');
cacheT2 = containers.Map('KeyType','double','ValueType','any');

for k = 1:numel(V)
    es = getCachedExp(cacheS2, ls(k), F.Zs);
    et = getCachedExp(cacheT2, lt(k), F.Zt);
    [es2, et2, ok] = applyMapsCell(es, et, modeS, toS, modeT, toT);
    if ~ok, continue; end
    kout2 = kout2 + 1;
    esStore(kout2,:) = es2(:).';
    etStore(kout2,:) = et2(:).';
    if kout2 == kout, break; end
end

% Build tensor indices from stored exponents using sorted bases
isH2 = zeros(kout,1);
itH2 = zeros(kout,1);
for kk = 1:kout
    subS = expVecToSubs(esStore(kk,:).', val2idxS);
    subT = expVecToSubs(etStore(kk,:).', val2idxT);
    isH2(kk) = sub2linTensor(subS, cellfun(@numel, ZsNew));
    itH2(kk) = sub2linTensor(subT, cellfun(@numel, ZtNew));
end

Inew = (brH-1)*dsN + isH2;
Jnew = (bcH-1)*dtN + itH2;

% Order triplets for faster sparse build
nRows = m*dsN;
lin = Inew + (Jnew-1)*nRows;
[~,p] = sort(lin);
Inew = Inew(p); Jnew = Jnew(p); vH = vH(p);

CH = sparse(Inew, Jnew, vH, m*dsN, n*dtN);
H  = quadPoly(CH, ZsNew, ZtNew, F.dim, F.ns, F.nt);

end

% ============================= Helpers =============================

function [vars, vals] = normalizePairs(vars, vals)
if ischar(vars) || isstring(vars), vars = {char(vars)}; end
if ischar(vals) || isstring(vals) || isnumeric(vals) || islogical(vals)
    vals = {vals};
end
if isstring(vars), vars = cellstr(vars); end
if isstring(vals), vals = cellstr(vals); end
if ~iscell(vars) || ~iscell(vals) || numel(vars) ~= numel(vals)
    error('subs: vars and vals must be same-length scalars/cells.');
end
for i = 1:numel(vars)
    vars{i} = char(string(vars{i}));
end
end

function [mp, k] = nameIndexMap(names)
% names: cellstr, sorted-unique
k = numel(names);
mp = containers.Map('KeyType','char','ValueType','double');
for i = 1:k
    mp(char(names{i})) = i;
end
end

function [modeS,toS,modeT,toT] = buildMaps(ks, kt, nsMap, ntMap, vars, vals)
% mode: 0 -> set to 0, 1 -> set to 1,
%       2 -> map to s(toS), 3 -> map to t(toS)  (for s-side vars)
%       2 -> map to s(toT), 3 -> map to t(toT)  (for t-side vars)

% defaults: identity
modeS = 2*ones(ks,1); toS = (1:ks).';   % s_i -> s_i
modeT = 3*ones(kt,1); toT = (1:kt).';   % t_i -> t_i (t/theta naming treated the same)

for k = 1:numel(vars)
    [side, idx] = parseVar(vars{k}, nsMap, ntMap);
    [rMode, rSide, rIdx] = parseRhs(vals{k}, nsMap, ntMap);

    if side == 's'
        if rMode <= 1
            modeS(idx) = rMode; toS(idx) = idx;
        elseif rSide == 's'
            modeS(idx) = 2; toS(idx) = rIdx;
        else
            modeS(idx) = 3; toS(idx) = rIdx;
        end
    else % t-side
        if rMode <= 1
            modeT(idx) = rMode; toT(idx) = idx;
        elseif rSide == 's'
            modeT(idx) = 2; toT(idx) = rIdx;
        else
            modeT(idx) = 3; toT(idx) = rIdx;
        end
    end
end
end

function [side, idx] = parseVar(v, nsMap, ntMap)
v = lower(strtrim(char(string(v))));
v = strrep(v, ' ', '');

% accept 's#' or any exact match in ns/nt
tok = regexp(v, '^(s|t|theta)(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    tag = tok{1};
    idxNum = str2double(tok{2});
    if strcmp(tag,'s')
        side = 's';
        idx = idxNum;
        return;
    else
        side = 't';
        idx = idxNum;
        return;
    end
end

% fallback: exact name match
if isKey(nsMap, v)
    side = 's'; idx = nsMap(v); return;
elseif isKey(ntMap, v)
    side = 't'; idx = ntMap(v); return;
else
    error('subs: bad variable "%s".', v);
end
end

function [rMode, rSide, rIdx] = parseRhs(x, nsMap, ntMap)
% rMode: 0 or 1 for constants; otherwise 2 (maps to s) or 3 (maps to t)
% rSide: 's' or 't' for mapping; rIdx target index
rSide = 's'; rIdx = 1;

if isnumeric(x) || islogical(x)
    rMode = (x~=0);
    return;
end

s = lower(strtrim(char(string(x))));
s = strrep(s, ' ', '');

if strcmp(s,'0'), rMode = 0; return; end
if strcmp(s,'1'), rMode = 1; return; end

tok = regexp(s, '^(s|t|theta)(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    tag = tok{1};
    rIdx = str2double(tok{2});
    if strcmp(tag,'s')
        rMode = 2; rSide = 's';
    else
        rMode = 3; rSide = 't';
    end
    return;
end

% exact match to existing variable names
if isKey(nsMap, s)
    rMode = 2; rSide = 's'; rIdx = nsMap(s);
elseif isKey(ntMap, s)
    rMode = 3; rSide = 't'; rIdx = ntMap(s);
else
    error('subs: RHS must be 0,1 or a variable name; got "%s".', s);
end
end

function [es2, et2, ok] = applyMapsCell(es, et, modeS, toS, modeT, toT)
% Apply substitution maps to exponent vectors (ks×1 and kt×1).
ks = numel(modeS);
kt = numel(modeT);

es2 = zeros(ks,1);
et2 = zeros(kt,1);
ok  = true;

% s-side variables
for i = 1:ks
    a = es(i);
    if a ~= 0
        ms = modeS(i);
        if ms == 0
            ok = false; return;              % s_i -> 0 kills term
        elseif ms == 1
            % s_i -> 1 : drop exponent
        elseif ms == 2
            es2(toS(i)) = es2(toS(i)) + a;   % s_i -> s_j
        else
            et2(toS(i)) = et2(toS(i)) + a;   % s_i -> t_j
        end
    end
end

% t-side variables
for i = 1:kt
    b = et(i);
    if b ~= 0
        mt = modeT(i);
        if mt == 0
            ok = false; return;              % t_i -> 0 kills term
        elseif mt == 1
            % t_i -> 1 : drop exponent
        elseif mt == 2
            es2(toT(i)) = es2(toT(i)) + b;   % t_i -> s_j
        else
            et2(toT(i)) = et2(toT(i)) + b;   % t_i -> t_j
        end
    end
end
end

function e = getCachedExp(cache, idx, Zcell)
key = double(idx);
if isKey(cache, key)
    e = cache(key);
    e = e(:);
else
    e = tensorExpAt(Zcell, idx);
    cache(key) = e;
end
end

function e = tensorExpAt(Zcell, idx)
% Decode tensor basis index to exponent vector; kron order, last varies fastest.
k = numel(Zcell);
e = zeros(k,1);
sz = cellfun(@numel, Zcell);
x = idx - 1;
for q = k:-1:1
    rq = sz(q);
    iq = mod(x, rq) + 1;
    e(q) = Zcell{q}(iq);
    x = floor(x / rq);
end
end

function [subs, Zsets, maps] = internExpVec(e, Zsets, maps)
% Insert each exponent component into per-variable set and return per-variable subscripts.
k = numel(e);
subs = zeros(k,1);
for i = 1:k
    val = double(e(i));
    mp = maps{i};
    if isKey(mp, val)
        subs(i) = mp(val);
    else
        Zsets{i} = [Zsets{i}; val]; %#ok<AGROW>
        idx = numel(Zsets{i});
        mp(val) = idx;
        maps{i} = mp;
        subs(i) = idx;
    end
end
end

function sz = sizesFromSets(Zsets)
sz = zeros(numel(Zsets),1);
for i = 1:numel(Zsets)
    sz(i) = max(1, numel(Zsets{i}));
end
end

function Zcell = finalizeSets(Zsets)
% Sort unique per variable; if empty, make [0].
k = numel(Zsets);
Zcell = cell(1,k);
for i = 1:k
    if isempty(Zsets{i})
        Zcell{i} = 0;
    else
        Zcell{i} = unique(Zsets{i}(:), 'sorted');
    end
end
end

function maps = valueToIndexMaps(Zcell)
k = numel(Zcell);
maps = cell(1,k);
for i = 1:k
    mp = containers.Map('KeyType','double','ValueType','double');
    z = Zcell{i}(:);
    for j = 1:numel(z)
        mp(double(z(j))) = j;
    end
    maps{i} = mp;
end
end

function subs = expVecToSubs(e, val2idx)
k = numel(e);
subs = zeros(k,1);
for i = 1:k
    subs(i) = val2idx{i}(double(e(i)));
end
end

function lin = sub2linTensor(subs, sizes)
% Convert per-variable subscripts to linear index in kron order (last varies fastest).
k = numel(sizes);
idx0 = 0;
stride = 1;
for i = k:-1:1
    idx0 = idx0 + (subs(i)-1)*stride;
    stride = stride * sizes(i);
end
lin = idx0 + 1;
end
