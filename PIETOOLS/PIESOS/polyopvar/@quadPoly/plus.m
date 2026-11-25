function H = plus(F, G)
%PLUS Add two quadPoly objects.
%
% INPUTS
%   F,G: quadPoly objects or 'double' matrices of the same matrix
%        dimension (not checked); a scalar is not expanded, except the
%        scalar 0, for which the other summand is returned
%
% OUTPUTS
%   H = F+G a quadPoly object, in the union of the monomial bases of F
%       and G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus
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

% Support addition with constant terms
if isequal(F,0)
    H = G;
    return
elseif isequal(G,0)
    H = F;
    return
end
if isa(F,'double')
    F = quadPoly(F);
end
if isa(G,'double')
    G = quadPoly(G);
end
% should add error checking (e.g. matching matrix dimensions)


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