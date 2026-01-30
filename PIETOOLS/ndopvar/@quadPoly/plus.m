function H = plus(F, G)
%PLUS Add two quadPoly objects under convention Z_old = P * Z_union.
%
% unionBasis signature:
%   [nU, ZU, PA, PB, mapA, mapB, unionIsA, unionIsB, dU] = unionBasis(nA, ZA, nB, ZB)
%
% We use the returned maps (old->union) for fast lifting.

% Fast path: identical representation
if isequal(F.ns,G.ns) && isequal(F.nt,G.nt) && isequal(F.Zs,G.Zs) && isequal(F.Zt,G.Zt)
    H = quadPoly(F.C + G.C, F.Zs, F.Zt, F.dim, F.ns, F.nt);
    return;
end

% --- s-side union ---
[nsU, ZsU, ~, ~, mapSF, mapSG, sIsF, sIsG, dsU] = unionBasis(F.ns, F.Zs, G.ns, G.Zs);

% --- t-side union ---
[ntU, ZtU, ~, ~, mapTF, mapTG, tIsF, tIsG, dtU] = unionBasis(F.nt, F.Zt, G.nt, G.Zt);

m = F.dim(1);
n = F.dim(2);

% If union equals F on both sides, lift only G
if sIsF && tIsF
    CG = liftIndexMaps(G.C, m, n, mapSG, mapTG, dsU, dtU);
    H  = quadPoly(F.C + CG, ZsU, ZtU, F.dim, nsU, ntU);
    return;
end

% If union equals G on both sides, lift only F
if sIsG && tIsG
    CF = liftIndexMaps(F.C, m, n, mapSF, mapTF, dsU, dtU);
    H  = quadPoly(CF + G.C, ZsU, ZtU, F.dim, nsU, ntU);
    return;
end

% Default: lift both
CF = liftIndex(F.C, m, n, mapSF, mapTF, dsU, dtU);
CG = liftIndex(G.C, m, n, mapSG, mapTG, dsU, dtU);

H  = quadPoly(CF + CG, ZsU, ZtU, F.dim, nsU, ntU);
end