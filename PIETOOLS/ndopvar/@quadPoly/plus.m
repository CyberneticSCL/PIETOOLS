function H = plus(F, G)
%plus Add two quadPoly objects
%
% Inputs:
%   F,G: quadPoly objects of same matrix dimension
%
% Outputs:
%   H = F+G a quadPoly object

% Fast path: if the monomials are same, just add coefficients
if isequal(F.ns,G.ns) && isequal(F.nt,G.nt) ...
        && isequal(F.Zs,G.Zs) && isequal(F.Zt,G.Zt)
    H = quadPoly(F.C + G.C, F.Zs, F.Zt, F.dim, F.ns, F.nt);
    return;
end

% --- left-monomial union ---
[nsU, ZsU, ~, ~, mapSF, mapSG, sIsF, sIsG, dsU] = unionBasis(F.ns, F.Zs, G.ns, G.Zs);

% --- right-monomial union ---
[ntU, ZtU, ~, ~, mapTF, mapTG, tIsF, tIsG, dtU] = unionBasis(F.nt, F.Zt, G.nt, G.Zt);

m = F.dim(1);
n = F.dim(2);

% If union equals F on both sides, lift only G
if sIsF && tIsF
    CG = liftIndex(G.C, m, n, mapSG, mapTG, dsU, dtU);
    H  = quadPoly(F.C + CG, ZsU, ZtU, F.dim, nsU, ntU);
    return;
end

% If union equals G on both sides, lift only F
if sIsG && tIsG
    CF = liftIndex(F.C, m, n, mapSF, mapTF, dsU, dtU);
    H  = quadPoly(CF + G.C, ZsU, ZtU, F.dim, nsU, ntU);
    return;
end

% Default: lift both
CF = liftIndex(F.C, m, n, mapSF, mapTF, dsU, dtU);
CG = liftIndex(G.C, m, n, mapSG, mapTG, dsU, dtU);

H  = quadPoly(CF + CG, ZsU, ZtU, F.dim, nsU, ntU);
end