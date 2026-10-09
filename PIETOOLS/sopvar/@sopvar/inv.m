function [Rinv,info] = inv(R,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [RINV,INFO] = INV(R[,OPTS]) the inverse of a square 1-D 'sopvar' block
%   (R x)(s) = R0(s) x(s) + int_a^s R1(s,t) x(t) dt + int_s^b R2(s,t) x(t) dt
% on L2^m[a,b], as a 'sopvar' with polynomial parameters of one degree d in
% s and in t. The block stores (sopvar document Sec. 4, 1-D)
%   R_i(s,t) = (I_m kron ZL(s)') C_i (I_m kron ZR(t)),   C_i = R.params{i+1},
% so each kernel is already in the separated form Lemma 16 below needs:
%   C_i = U_i V_i  (SVD, numerical rank r_i)  ==>  R_i = F_i(s) G_i(t),
%   F_i(s) = (I kron ZL(s)') U_i   (m x r_i),   G_i(t) = V_i (I kron ZR(t)).
%
% METHOD (Shivakumar, Das, Peet, arXiv 2208.13104, Sec. VII, Lemma 16 and
% Cor. 17, Gohberg-Krein). With H_i = R0^{-1} R_i = -F_i^h G_i, F_i^h =
% -R0^{-1} F_i, set C(s) = [F_1^h F_2^h] (m x n), B(s) = [G_1; -G_2] (n x m),
% n = r_1 + r_2, and solve on [a,b]
%   U' = B C U, U(a) = I;      V' = -V B C, V(a) = I   (V = U^{-1}).
% Partition U(b) with U22 of size r_2; R is invertible iff U22 is, and
%   Pg = [0 0; U22^{-1} U21 I],
%   R^{-1} = {R0^{-1},  C(s) U(s) (I - Pg) V(t) B(t) R0(t)^{-1},
%                      -C(s) U(s) Pg V(t) B(t) R0(t)^{-1}}.
% GK_GRID (sopvar/misc) propagates U, V by RK4 on OPTS.N uniform nodes with
% midpoints and returns M(s) = C(s)U(s), N(t) = V(t)B(t)R0(t)^{-1} and
% R0^{-1} at the nodes. The three parameters are evaluated on the node grid
% (the kernels on the WHOLE square from the formulas, not only on their
% triangle, because the formula is the smooth continuation and a
% polynomial fits it; fitting a kernel that is cut at the diagonal does not
% converge) and fitted by least squares in the Chebyshev basis on [a,b]
% (CHEBFIT_MONO), one degree d for the three, raised from OPTS.deg0 in
% steps of 2 until the relative RMS residual of all three fits is below
% OPTS.tol or d = OPTS.degmax; the monomial coefficients are placed on
% ZL = ZR = (0:d)', the multiplier in the ZR-degree-0 columns, the canonical
% multiplier form of the class. R0^{-1} is rational, so the inverse is only
% approximated by a polynomial block; INFO reports the fit residuals. This
% is the route of the stock 'inv_opvar_2' (SS, 02/2026), written against
% the block storage: no polynomial objects, no monomial bookkeeping, the
% kernel split is one SVD per kernel.
%
% INPUTS
% - R:      'sopvar', L2^m[s] -> L2^m[s] with one spatial variable, square;
% - opts:   (optional) struct, fields
%     N        nodes of the RK4 grid, default 101 (the midpoints are added);
%     tol      target relative RMS fit residual, default 1e-8;
%     deg0     first fit degree, default 4;
%     degmax   last fit degree tried, default 16;
%     deg      a fixed fit degree (overrides deg0, degmax, tol);
%     tol_rank relative SVD cut for the kernel rank, default 1e-13.
% OUTPUTS
% - Rinv:   'sopvar', the inverse block on ZL = ZR = (0:d)';
% - info:   struct: N, h, r1, r2 (kernel ranks), n, rcondU22, condR0 (max over
%           the grid), d, relrms (1 x 3, the fit residuals of R0^{-1}, R1^h,
%           R2^h), and note.
%
% NOTES
% Cost: O(N n^3) for the RK4 sweep, O(N^2 m^2 n) for the kernel values,
% O(N^2 (d+1)^2 m^2) per fit degree. No decision variables are involved; the
% block is fixed, so nothing here scales with q.
% A singular R0(s) at a grid point, or a near-singular U22(b), is reported
% through INFO and a warning; the formulas then give a large residual rather
% than an error.
%
% See also COPVAR/INV, GK_GRID, CHEBFIT_MONO, GETCONTROLLER_DIRECT_SOP,
% INV_OPVAR_2, INV_OPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - inv (sopvar)
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
% Initial coding MMP, 10/09/2026. Container counterpart of @opvar/inv for
%                the synthesis reconstruction K = Z P^{-1}, L = P^{-1} Z on
%                the container path (getController_sop, getObserver_sop).
% MMP, 10/09/2026: The grid part (kernel factors, RK4, Pg, M, N) moved to
%                sopvar/misc/gk_grid.m and the Chebyshev fit to
%                sopvar/misc/chebfit_mono.m, both shared with
%                getController_direct_sop, which builds the gains from the
%                same grid data without forming the inverse. Same output
%                (Test_copvar_inv passes unchanged).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    opts = struct();    end
o = struct('N',101,'tol',1e-8,'deg0',4,'degmax',16,'deg',[],'tol_rank',1e-13);
fn = fieldnames(opts);
for k = 1:numel(fn),    o.(fn{k}) = opts.(fn{k});   end
if ~isempty(o.deg),     o.deg0 = o.deg;     o.degmax = o.deg;   end

% ---- the grid data (validates the block)
G = gk_grid(R,struct('N',o.N,'tol_rank',o.tol_rank));
m = G.m;    N = G.N;    n = G.n;    a = G.a;    b = G.b;    sg = G.sg;
Pg = G.Pg;  ImP = eye(n) - Pg;

% ---- parameter values on the node grid: S0_i = R0^{-1}(s_i),
%      S1_ij = M_i (I-Pg) N_j, S2_ij = -M_i Pg N_j; whole square, see header
S1 = zeros(m,m,N,N);    S2 = zeros(m,m,N,N);
if n>0
    for i = 1:N
        G1 = pagemtimes(G.M(:,:,i)*ImP,G.Nn);           % m x m x N over j
        G2 = pagemtimes(-G.M(:,:,i)*Pg,G.Nn);
        S1(:,:,i,:) = reshape(G1,m,m,1,N);
        S2(:,:,i,:) = reshape(G2,m,m,1,N);
    end
end

% ---- fit: one degree d for all three, raised until every residual meets tol
Y0 = reshape(permute(G.R0inv,[3 1 2]),N,m*m);           % N x m^2, column-major entries
Y1 = reshape(permute(S1,[3 4 1 2]),N*N,m*m);            % rows (i,j), i fastest
Y2 = reshape(permute(S2,[3 4 1 2]),N*N,m*m);
d = o.deg0;
while true
    [C0mono,e0] = chebfit_mono(sg,Y0,[a b],d,d,0);
    if n>0
        [C1mono,e1] = chebfit_mono(sg,Y1,[a b],d,d,0,sg);
        [C2mono,e2] = chebfit_mono(sg,Y2,[a b],d,d,0,sg);
    else
        C1mono = zeros((d+1)^2,m*m);    C2mono = C1mono;    e1 = 0;    e2 = 0;
    end
    if max([e0 e1 e2])<=o.tol || d>=o.degmax,   break,  end
    d = min(d+2,o.degmax);
end

% ---- place into the block storage on ZL = ZR = (0:d)'
nZ = d + 1;
C0h = zeros(m*nZ,m*nZ);     C1h = C0h;     C2h = C0h;
for i = 1:m
    for j = 1:m
        col = i + (j-1)*m;                              % column-major entry (i,j)
        C0h((i-1)*nZ+(1:nZ),(j-1)*nZ+1) = C0mono(:,col);        % ZR degree 0 column: canonical multiplier form
        C1h((i-1)*nZ+(1:nZ),(j-1)*nZ+(1:nZ)) = reshape(C1mono(:,col),nZ,nZ).';   % (r,u) -> row r, col u
        C2h((i-1)*nZ+(1:nZ),(j-1)*nZ+(1:nZ)) = reshape(C2mono(:,col),nZ,nZ).';
    end
end
Zh = (0:d)';
Rinv = sopvar({sparse(C0h);sparse(C1h);sparse(C2h)},R.vars,{Zh},{Zh},R.dom,[m m]);

info = struct('N',N,'h',G.h,'r1',G.r1,'r2',G.r2,'n',n,'rcondU22',G.rcondU22,'condR0',G.condR0,...
              'd',d,'relrms',[e0 e1 e2],'note','Lemma 16 / Cor. 17 of arXiv 2208.13104, RK4 + Chebyshev LS fit');
if max(info.relrms)>o.tol
    info.note = sprintf('%s; fit tolerance %.1e NOT met at degmax %d (relrms %.2e)',info.note,o.tol,o.degmax,max(info.relrms));
end
end
