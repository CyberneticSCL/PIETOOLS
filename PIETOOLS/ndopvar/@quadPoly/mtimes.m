function H = mtimes(A,B)
%MTIMES Multiplication for quadPoly.
%
% Supports:
%   quadPoly * scalar double
%   scalar double * quadPoly
%   constant double matrix * quadPoly
%   quadPoly * constant double matrix
%   quadPoly * quadPoly

% ---------- scalar double multiply ----------
if isa(A,'quadPoly') && isa(B,'double') && isscalar(B)
    H = quadPoly(A.C * B, A.Zs, A.Zt, A.dim, A.ns, A.nt);
    return;
elseif isa(A,'double') && isscalar(A) && isa(B,'quadPoly')
    H = quadPoly(A * B.C, B.Zs, B.Zt, B.dim, B.ns, B.nt);
    return;
end

% ---------- constant matrix multiply ----------
if isa(A,'double') && isa(B,'quadPoly')
    H = leftConstTimesPoly(A, B);
    return;
elseif isa(A,'quadPoly') && isa(B,'double')
    H = rightConstTimesPoly(A, B);
    return;
end

% ---------- quadPoly * quadPoly ----------
if ~(isa(A,'quadPoly') && isa(B,'quadPoly'))
    error('quadPoly:mtimes','Unsupported operands for mtimes.');
end

F = A; 
G = B;

% scalar-poly scaling (matrix * scalar poly)
if all(G.dim == [1 1])
    H = scalarPolyScale(F, G); % F .* G (entrywise poly multiply)
    return;
elseif all(F.dim == [1 1])
    H = scalarPolyScale(G, F); % (scalar) .* G
    return;
end

% dimension checks
m = F.dim(1); n = F.dim(2); p = G.dim(2);
if n ~= G.dim(1)
    error('quadPoly:mtimes','Inner dimensions must match.');
end

nsF = size(F.Zs,2); nsG = size(G.Zs,2);
ntF = size(F.Zt,2); ntG = size(G.Zt,2);
if nsF ~= nsG || ntF ~= ntG
    error('quadPoly:mtimes','Left/right variable counts must match for multiplication.');
end

% decode sparse coefficient matrices
dsF = size(F.Zs,1); dtF = size(F.Zt,1);
dsG = size(G.Zs,1); dtG = size(G.Zt,1);

[iF,aF,jF,bF,vF] = decodeCoeff(F.C, dsF, dtF);
if isempty(vF)
    H = zeroQuadPoly([m p], nsF, ntF, F.ns, F.nt);
    return;
end

[jG,cG,kG,dG,vG] = decodeCoeff(G.C, dsG, dtG); % note: row-block is "j"
if isempty(vG)
    H = zeroQuadPoly([m p], nsF, ntF, F.ns, F.nt);
    return;
end

% ---------- group by shared index j ----------
[jF_sorted, ordF] = sort(jF);
iF = iF(ordF); aF = aF(ordF); bF = bF(ordF); vF = vF(ordF);

[jG_sorted, ordG] = sort(jG);
cG = cG(ordG); kG = kG(ordG); dG = dG(ordG); vG = vG(ordG);

ptrF = groupPtrFromSorted(jF_sorted, n);
ptrG = groupPtrFromSorted(jG_sorted, n);

% ---------- build output bases + coefficient triplets ----------
ZsH = zeros(0, nsF);
ZtH = zeros(0, ntF);

expMapS = containers.Map('KeyType','char','ValueType','double');
expMapT = containers.Map('KeyType','char','ValueType','double');

% cache basis indices for (a,c) pairs and (b,d) pairs
pairS = containers.Map('KeyType','uint64','ValueType','double');
pairT = containers.Map('KeyType','uint64','ValueType','double');

% triplet blocks (collected by chunks to avoid constant realloc)
Iblk = {}; Jblk = {}; Vblk = {};
tcount = 0;

dsKeyStride = uint64(dsF);
dtKeyStride = uint64(dtF);

for j = 1:n
    f0 = ptrF(j);    f1 = ptrF(j+1)-1;
    g0 = ptrG(j);    g1 = ptrG(j+1)-1;
    if f0 > f1 || g0 > g1
        continue;
    end

    % F terms for this j: (i,a,b,v)
    iFj = iF(f0:f1);
    aFj = aF(f0:f1);
    bFj = bF(f0:f1);
    vFj = vF(f0:f1);

    % G terms for this j: (c,k,d,v)
    cGj = cG(g0:g1);
    kGj = kG(g0:g1);
    dGj = dG(g0:g1);
    vGj = vG(g0:g1);

    ng = numel(vGj);

    for u = 1:numel(vFj)
        ii = iFj(u);
        aa = aFj(u);
        bb = bFj(u);
        vv = vFj(u);

        % compute/lookup all needed s-basis indices for (aa, cGj(:))
        isVec = zeros(ng,1);
        keyS = uint64(aa) + (uint64(cGj(:))-1)*dsKeyStride;
        for t = 1:ng
            ks = keyS(t);
            if isKey(pairS, ks)
                isVec(t) = pairS(ks);
            else
                e = F.Zs(aa,:) + G.Zs(cGj(t),:);
                [idx, ZsH] = getOrAddExp(expMapS, e, ZsH);
                pairS(ks) = idx;
                isVec(t) = idx;
            end
        end

        % compute/lookup all needed t-basis indices for (bb, dGj(:))
        itVec = zeros(ng,1);
        keyT = uint64(bb) + (uint64(dGj(:))-1)*dtKeyStride;
        for t = 1:ng
            kt = keyT(t);
            if isKey(pairT, kt)
                itVec(t) = pairT(kt);
            else
                e = F.Zt(bb,:) + G.Zt(dGj(t),:);
                [idx, ZtH] = getOrAddExp(expMapT, e, ZtH);
                pairT(kt) = idx;
                itVec(t) = idx;
            end
        end

        % append triplets for this (F-term) x (all G-terms)
        tcount = tcount + 1;

        dsHtmp = max(1, size(ZsH,1));
        dtHtmp = max(1, size(ZtH,1));

        Iblk{tcount,1} = (ii-1)*dsHtmp + isVec;
        Jblk{tcount,1} = (kGj(:)-1)*dtHtmp + itVec;
        Vblk{tcount,1} = vv .* vGj(:);
    end
end

% finalize basis (ensure at least one monomial)
if isempty(ZsH), ZsH = zeros(1,nsF); end
if isempty(ZtH), ZtH = zeros(1,ntF); end

dsH = size(ZsH,1);
dtH = size(ZtH,1);

if tcount == 0
    CH = sparse(m*dsH, p*dtH);
else
    I = vertcat(Iblk{:});
    J = vertcat(Jblk{:});
    V = vertcat(Vblk{:});
    CH = sparse(I, J, V, m*dsH, p*dtH); % sums duplicates automatically
end

H = quadPoly(CH, ZsH, ZtH, [m p], F.ns, F.nt);

end

% ===================== helpers =====================

function H = leftConstTimesPoly(L, F)
% (double) * (quadPoly)
ds = size(F.Zs,1);
Cnew = kron(sparse(L), speye(ds)) * F.C;
H = quadPoly(Cnew, F.Zs, F.Zt, [size(L,1) F.dim(2)], F.ns, F.nt);
end

function H = rightConstTimesPoly(F, R)
% (quadPoly) * (double)
dt = size(F.Zt,1);
Cnew = F.C * kron(sparse(R), speye(dt));
H = quadPoly(Cnew, F.Zs, F.Zt, [F.dim(1) size(R,2)], F.ns, F.nt);
end

function H = zeroQuadPoly(dim, ns, nt, namesS, namesT)
% canonical zero polynomial with 1 monomial in each basis
Zs = zeros(1, ns);
Zt = zeros(1, nt);
C  = sparse(dim(1), dim(2)); % since ds=dt=1
H  = quadPoly(C, Zs, Zt, dim, namesS, namesT);
end

function [i,a,j,b,v] = decodeCoeff(C, ds, dt)
% Decode sparse C into block indices:
% row = (i-1)*ds + a, col = (j-1)*dt + b
[I,J,v] = find(C);
if isempty(v)
    i=[]; a=[]; j=[]; b=[]; return;
end
i = floor((I-1)/ds) + 1;
a = I - (i-1)*ds;
j = floor((J-1)/dt) + 1;
b = J - (j-1)*dt;
end

function ptr = groupPtrFromSorted(js_sorted, n)
% js_sorted is sorted and in 1..n
ptr = ones(n+1,1);
if isempty(js_sorted)
    ptr(:) = 1;
    return;
end
counts = accumarray(js_sorted, 1, [n 1]);
ptr(1) = 1;
for k = 1:n
    ptr(k+1) = ptr(k) + counts(k);
end
end

function [idx, Z] = getOrAddExp(mp, e, Z)
% Store exponent row e in Z (if new) and return its 1-based row index.
% Key as comma-separated ints (simple + robust).
key = sprintf('%d,', e);
if isKey(mp, key)
    idx = mp(key);
else
    Z(end+1,:) = e;
    idx = size(Z,1);
    mp(key) = idx;
end
end

function H = scalarPolyScale(F, S)
% Entrywise polynomial multiplication of matrix poly F by scalar poly S (1x1).
% Result has same dim as F.
m = F.dim(1); n = F.dim(2);

ns = size(F.Zs,2); nt = size(F.Zt,2);
if size(S.Zs,2) ~= ns || size(S.Zt,2) ~= nt
    error('quadPoly:mtimes','Scalar poly variable counts must match.');
end

dsF = size(F.Zs,1); dtF = size(F.Zt,1);
dsS = size(S.Zs,1); dtS = size(S.Zt,1);

[iF,aF,jF,bF,vF] = decodeCoeff(F.C, dsF, dtF);
if isempty(vF)
    H = zeroQuadPoly([m n], ns, nt, F.ns, F.nt);
    return;
end

% S is 1x1, so indices are directly (aS,bS) in its dsS x dtS C matrix
[aS,bS,vS] = find(S.C);
if isempty(vS)
    H = zeroQuadPoly([m n], ns, nt, F.ns, F.nt);
    return;
end

ZsH = zeros(0, ns);
ZtH = zeros(0, nt);
expMapS = containers.Map('KeyType','char','ValueType','double');
expMapT = containers.Map('KeyType','char','ValueType','double');

pairS = containers.Map('KeyType','uint64','ValueType','double');
pairT = containers.Map('KeyType','uint64','ValueType','double');

dsKeyStride = uint64(dsF);
dtKeyStride = uint64(dtF);

Iblk = {}; Jblk = {}; Vblk = {};
tcount = 0;

for u = 1:numel(vF)
    ii = iF(u);
    jj = jF(u);
    aa = aF(u);
    bb = bF(u);
    vv = vF(u);

    nsTerms = numel(vS);
    isVec = zeros(nsTerms,1);
    itVec = zeros(nsTerms,1);

    keyS = uint64(aa) + (uint64(aS)-1)*dsKeyStride;
    keyT = uint64(bb) + (uint64(bS)-1)*dtKeyStride;

    for t = 1:nsTerms
        ks = keyS(t);
        if isKey(pairS, ks)
            isVec(t) = pairS(ks);
        else
            e = F.Zs(aa,:) + S.Zs(aS(t),:);
            [idx, ZsH] = getOrAddExp(expMapS, e, ZsH);
            pairS(ks) = idx;
            isVec(t) = idx;
        end

        kt = keyT(t);
        if isKey(pairT, kt)
            itVec(t) = pairT(kt);
        else
            e = F.Zt(bb,:) + S.Zt(bS(t),:);
            [idx, ZtH] = getOrAddExp(expMapT, e, ZtH);
            pairT(kt) = idx;
            itVec(t) = idx;
        end
    end

    % append triplets
    tcount = tcount + 1;

    dsHtmp = max(1, size(ZsH,1));
    dtHtmp = max(1, size(ZtH,1));

    Iblk{tcount,1} = (ii-1)*dsHtmp + isVec;
    Jblk{tcount,1} = (jj-1)*dtHtmp + itVec;
    Vblk{tcount,1} = vv .* vS(:);
end

if isempty(ZsH), ZsH = zeros(1,ns); end
if isempty(ZtH), ZtH = zeros(1,nt); end

dsH = size(ZsH,1);
dtH = size(ZtH,1);

I = vertcat(Iblk{:});
J = vertcat(Jblk{:});
V = vertcat(Vblk{:});

CH = sparse(I, J, V, m*dsH, n*dtH);
H  = quadPoly(CH, ZsH, ZtH, [m n], F.ns, F.nt);
end
