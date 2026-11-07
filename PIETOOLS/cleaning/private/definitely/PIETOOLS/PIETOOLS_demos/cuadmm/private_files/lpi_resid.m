function r = lpi_resid(prog,M,Atf,bf,D)
% lpi_resid -- the full residual panel for a solved stability LPI.
%
% Deliberately PLURAL.  There is no single right number here, and the point is
% to see how the measures relate, not to pick one:
%
%  ALGEBRAIC (two conventions, both reported)
%    rel_b  = ||A'x-b|| / ||b||        what PIETOOLS/SeDuMi report
%    rel_cu = ||A'x-b|| / (1+||b||)    what cuADMM's stopping rule uses
%                                      (solver.cu:218; equals rel_b/2 once
%                                       sossolve has normalised ||b||=1)
%
%  CONE
%    psd_min    smallest eigenvalue over all Gram blocks
%    psd_relmin the same, scaled by each block's largest eigenvalue
%
%  OPERATOR, from the Galerkin matrix of the LPI residual R = Dop + Deop
%    op_ind   ||R||_2   induced L2->L2, relative to ||Dop||_2
%    op_hs    ||R||_F   Hilbert-Schmidt, relative to ||Dop||_F
%    op_coef  coefficient-space Frobenius over parameter cells (what opcheck
%             does) -- a different weighting of the same information
%
%  POSITIVITY AS AN OPERATOR (not merely as a Gram)
%    pop_mineig  smallest eigenvalue of the Galerkin matrix of Pop.  A PSD
%                Gram does not by itself make the OPERATOR positive, so this
%                is checked separately.

if nargin<5 || isempty(D), D = 16; end

x  = prog.solinfo.RRx(:);
res = Atf'*x - bf;
nb  = norm(full(bf));
r.alg_abs = norm(full(res));
r.rel_b   = r.alg_abs / nb;
r.rel_cu  = r.alg_abs / (1 + nb);
r.normb   = nb;

% ---- cone: every Gram block
% Normalise the PSD margin by the LARGEST eigenvalue ACROSS ALL BLOCKS, not by
% each block's own maximum.  MEASURED FAILURE (2026-09-23): a numerically zero
% block has eigenvalues that are pure rounding noise, so its own-max ratio is
% ~-1 regardless of the solution, and the per-block version reported -1.000 at
% every eppos2 -- a property of the normalisation, not of the solver.
S = sdpshape(prog);
r.psd_min = inf;  r.psd_max = -inf;  r.psd_blocks = zeros(numel(S.Ks),2);
off = S.Kf;
for k = 1:numel(S.Ks)
    N  = S.Ks(k);
    Xk = reshape(x(off+(1:N^2)), N, N);
    Xk = (Xk+Xk')/2;
    ev = eig(Xk);
    r.psd_blocks(k,:) = [min(ev) max(ev)];
    r.psd_min = min(r.psd_min, min(ev));
    r.psd_max = max(r.psd_max, max(ev));
    off = off + N^2;
end
r.psd_relmin = r.psd_min / max(r.psd_max, eps);

% ---- operator level
Dsol  = getsol_lpivar(prog, M.Dop);
Desol = getsol_lpivar(prog, M.Deop);
Rres  = Dsol + Desol;                     % must be the zero operator

GR = opnorm_pi(Rres,D);
GD = opnorm_pi(Dsol,D);
r.op_ind  = norm(GR)          / max(norm(GD),eps);
r.op_hs   = norm(GR,'fro')    / max(norm(GD,'fro'),eps);
r.op_coef = coefnorm(Rres)    / max(coefnorm(Dsol),eps);
r.D       = D;

% ABSOLUTE scale of the thing being normalised against, and a triviality flag.
% MEASURED TRAP (2026-09-23): on an infeasible program SeDuMi still returns a
% solution vector, and it is X=0.  Then Dsol=0, so op_ind and op_hs are 0/0 ->
% 0 and psd_min is exactly 0: the operator residual and the cone both look
% PERFECT.  Only the algebraic residual exposes it (rel_b = 1 exactly, since
% A'*0-b = -b).  A gate on the operator residual alone accepts X=0.
r.normDop = norm(GD);                 % absolute; op_ind is meaningless if ~0
r.normx   = norm(x);
% Trivial iff the solution vector is (numerically) zero, equivalently the
% negativity operator vanished, equivalently the residual equals -b.
r.trivial = (r.normDop <= 1e-12) || (abs(r.rel_b - 1) <= 1e-6) || (r.normx <= 1e-12);

% ---- is the Lyapunov operator positive AS AN OPERATOR?
Psol = getsol_lpivar(prog, M.Pop);
GP   = opnorm_pi(Psol,D);
GPs  = (GP+GP')/2;
r.pop_mineig = min(eig(GPs));
r.pop_norm   = norm(GP);

% Operator residual against a denominator that CANNOT COLLAPSE.  ||Dop|| is
% the natural scale but it is itself linear in the decision variables, so at
% X=0 it is zero and op_ind reads 0/0 -> 0, i.e. a PERFECT residual for the
% trivial point (measured).  Pop carries the CONSTANT eppos term, so ||Pop||
% is bounded away from zero whatever the solver returns.  Gate on this, not on
% op_ind.  (Convention adopted from the low-rank package's opcheck, which
% records the same failure: "rel is INVALID wherever the candidate drives
% Dop -> 0 ... Accept on relP, never on rel alone".)
r.op_indP = norm(GR) / max(r.pop_norm,eps);
r.op_hsP  = norm(GR,'fro') / max(norm(GP,'fro'),eps);

% THE OTHER OPERATOR EQUALITY.  PIE2PDEstability imposes TWO: Top'*Qop == Pop
% and Dop + Deop == 0.  Dop = A'Q + Q'A is entirely linear in the decision
% variable Q -- the eppos constant sits in Pop, which appears only in the
% FIRST equality -- so at X=0 the NEGATIVITY residual is exactly zero,
% numerator and denominator alike, and NO choice of denominator on that
% constraint can reject the trivial point (measured: op_ind and op_indP both
% 0.00e+00 there).  The first equality does reject it: T'*0 - Pop = -eppos*T'T.
% LESSON: an operator-level gate must cover EVERY operator equality the
% program imposes, not just the one the analysis is "about".
if isfield(M,'Qop') && ~isempty(M.Qop)
    Qsol  = getsol_lpivar(prog, M.Qop);
    GEQ   = opnorm_pi(M.Top'*Qsol - Psol, D);
    r.op_eqQ = norm(GEQ) / max(r.pop_norm,eps);
else
    r.op_eqQ = NaN;
end
end


function v = coefnorm(Op)
% Frobenius norm over every opvar parameter cell, in COEFFICIENT space.
% This is the measure opcheck uses; it weights the cells by however the
% operator happens to be parameterised, which is why it is reported next to
% the two basis-free norms rather than instead of them.
v = 0;
f = {'P','Q1','Q2'};
for i = 1:numel(f)
    v = v + cellsq(Op.(f{i}));
end
g = {'R0','R1','R2'};
for i = 1:numel(g)
    v = v + cellsq(Op.R.(g{i}));
end
v = sqrt(v);
end

function t = cellsq(A)
t = 0;
if isempty(A), return; end
if isa(A,'double'), t = sum(A(:).^2); return; end
c = A.coefficient;
if ~isempty(c), t = full(sum(sum(c.^2))); end
end
