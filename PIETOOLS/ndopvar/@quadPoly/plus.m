function H = plus(F, G)
% Fast path: identical bases -> just add sparse matrices
if isequal(F.Zs, G.Zs) && isequal(F.Zt, G.Zt)
    H = quadPoly(F.C + G.C, F.Zs, F.Zt, F.dim, F.ns, F.nt);
    return;
end

% Union bases + maps
[ZsU, mapSF, mapSG] = unionBasis(F.Zs, G.Zs);
[ZtU, mapTF, mapTG] = unionBasis(F.Zt, G.Zt);

m = F.dim(1); n = F.dim(2);
dsF = size(F.Zs,1); dtF = size(F.Zt,1);
dsG = size(G.Zs,1); dtG = size(G.Zt,1);
dsU = size(ZsU,1); dtU = size(ZtU,1);

% If the union is exactly F's basis, only lift G once
if dsU==dsF && dtU==dtF && isequal(mapSF,(1:dsF)') && isequal(mapTF,(1:dtF)')
    CG = lift(G.C, m, n, dsG, dtG, dsU, dtU, mapSG, mapTG);
    H  = quadPoly(F.C + CG, ZsU, ZtU, F.dim, F.ns, F.nt);
    return;
end

% If the union is exactly G's basis, only lift F once
if dsU==dsG && dtU==dtG && isequal(mapSG,(1:dsG)') && isequal(mapTG,(1:dtG)')
    CF = lift(F.C, m, n, dsF, dtF, dsU, dtU, mapSF, mapTF);
    H  = quadPoly(CF + G.C, ZsU, ZtU, F.dim, F.ns, F.nt);
    return;
end

% Default: lift both and add (typically fastest/most stable in MATLAB)
CF = lift(F.C, m, n, dsF, dtF, dsU, dtU, mapSF, mapTF);
CG = lift(G.C, m, n, dsG, dtG, dsU, dtU, mapSG, mapTG);
H  = quadPoly(CF + CG, ZsU, ZtU, F.dim, F.ns, F.nt);
end

