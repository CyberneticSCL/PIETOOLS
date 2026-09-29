function [keep,info] = heatNd_lindep(At,rtol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [KEEP,INFO] = HEATND_LINDEP(AT,RTOL) a maximal linearly independent set
% of the equality rows of the SDP At'x = b (the columns of At), by one
% Q-less sparse QR with column pivoting (SuiteSparseQR, via 'qr'): a column
% whose remaining norm falls below the tolerance is dead, so the rows of R
% with |R(j,j)| > RTOL * max|diag R| index an independent set.
%
% WHY. The heatNd equality systems are rank-deficient (measured 2-D: rank
% 990 of 1242 rows) and MOSEK's presolve removes none ("Lin. dep. - primal
% deps.: 0"), so its Schur complement is singular. HEATND_SOLVE can solve
% on the rows KEEP only; the verdict is still computed on the FULL system,
% so a wrongly dropped row can only turn a verdict uncertain, never wrong.
% The dependency is structural (kernel-coefficient identities), so a set
% from one generic k is reused across a bisection; HEATND_SOLVE rechecks.
%
% INPUT  At: nx x m sparse; rtol: default 1e-10.
% OUTPUT keep: sorted row indices; info: rank, m, t (s), rtol.
% Cost: one sparse QR of At (fill ~ the Cholesky of At'At, the same order
% as ONE interior-point factorization); nothing dense.
%
% Initial coding MMP, 09/27/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(rtol),   rtol = 1e-10;   end
t0 = tic;
% Scale columns to unit norm first, so the tolerance is relative per row.
cn = full(sqrt(sum(At.^2,1)));  cn(cn==0) = 1;
As = At*spdiags(1./cn(:),0,numel(cn),numel(cn));
[~,Rq,p] = qr(As,sparse(size(As,1),0),0);       % Q-less, 'vector' permutation
dg = abs(full(diag(Rq)));
keep = sort(p(dg > rtol*max(dg)));
keep = keep(:);
info = struct('rank',numel(keep),'m',size(At,2),'t',toc(t0),'rtol',rtol);
end
