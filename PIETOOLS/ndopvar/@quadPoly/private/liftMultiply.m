function Cup = liftMultiply(C, m, n, Ps, Pt)
%LIFTMULTIPLY Lift coefficient matrix via reshape-based multiplies (no kron).
%
% Convention (selection/basis maps):
%   Z_old = P * Z_union
% where P is a row-pick (selection) matrix of size d_old × d_union.
%
% For F(s,t) = (I_m ⊗ Zs_old^T) * C * (I_n ⊗ Zt_old),
% with Zs_old = Ps * Zs_union  and  Zt_old = Pt * Zt_union,
% the lifted coefficient matrix in the union basis is:
%   Cup = (I_m ⊗ Ps') * C * (I_n ⊗ Pt)
%
% Sizes:
%   C  : (m*dsOld) × (n*dtOld)
%   Ps : dsOld × dsU      (s-side selection; old = Ps * union)
%   Pt : dtOld × dtU      (t-side selection; old = Pt * union)
%
% Output:
%   Cup: (m*dsU) × (n*dtU)

C  = sparse(C);
Ps = sparse(Ps);
Pt = sparse(Pt);

[dsOld, dsU] = size(Ps);
[dtOld, dtU] = size(Pt);

% ---- Left multiply by (I_m ⊗ Ps') ----
% Interpret C as m row-blocks of size dsOld, across all columns.
X  = reshape(C, dsOld, []);          % dsOld × (m*n*dtOld)
Y  = Ps.' * X;                       % dsU   × (m*n*dtOld)
C1 = reshape(Y, m*dsU, n*dtOld);     % (m*dsU) × (n*dtOld)

% ---- Right multiply by (I_n ⊗ Pt) ----
% Right multiplication mixes dt-blocks within each of the n column-blocks.
% Work on transpose: right-multiply by (I_n ⊗ Pt) <=> left-multiply C1' by (I_n ⊗ Pt')
Xt  = reshape(C1.', dtOld, []);      % dtOld × (n*m*dsU)
Yt  = Pt.' * Xt;                     % dtU   × (n*m*dsU)
Cup = reshape(Yt, n*dtU, m*dsU).';   % (m*dsU) × (n*dtU)

end
