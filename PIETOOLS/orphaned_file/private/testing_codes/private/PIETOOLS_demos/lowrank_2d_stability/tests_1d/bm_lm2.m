function [w,f,it] = bm_lm2(P,rv,w0,maxit,tol)
% bm_lm2 -- Levenberg-Marquardt for Burer-Monteiro, choosing the CHEAPER normal
% equations.  Drop-in replacement for bm_lm; identical maths, different pivot.
%
% WHY.  bm_lm always solves in the DUAL m-space:  d = -J'(JJ' + lam I)^{-1} F.
% That is the right choice in 1D, where the comment "the system is
% underdetermined (m << numel(w))" holds.  In 2D AT LOW RANK it is inverted:
%   light n=1, rv = [1 1]:  nv = sum(Ns.*rv) = 432  vs  mres = rank(A) = 2689
% so the system is OVERdetermined and the dual form factorises a 2689 x 2689
% matrix (~7 GFlop of Cholesky, inside a 40-try damping loop) to take a step in
% a 432-dimensional space.  The primal form
%   d = -(J'J + lam I)^{-1} J'F
% factorises 432 x 432 instead -- the same Gauss-Newton step, ~240x less
% factorisation work.  Low rank is exactly the regime under test here, so the
% dual form would have made the cross-check unaffordable rather than wrong.
%
% Both branches are kept and selected on the measured shapes, so this stays
% correct if a later caller runs at high rank.
if nargin<4||isempty(maxit), maxit=200; end
if nargin<5||isempty(tol),   tol=1e-13;  end
w = w0;  lam = 1e-8;
[F,J] = bm_resid(w,P,rv);  f = norm(F);
for it = 1:maxit
    nv = size(J,2);   mr = size(J,1);
    primal = nv <= mr;
    if primal
        A = full(J'*J);   g = J'*F;
    else
        A = full(J*J');
    end
    sc = max(max(abs(diag(A))),eps);
    ok = false;
    for tries = 1:40
        M = A + (lam*sc)*eye(size(A,1));
        [L,p] = chol(M,'lower');
        if primal
            if p==0, d = -(L'\(L\g)); else, d = -(M\g); end
        else
            if p==0, y = L'\(L\F); else, y = M\F; end
            d = -(J'*y);
        end
        wn = w + d;
        Fn = bm_resid(wn,P,rv);   fn = norm(Fn);
        if isfinite(fn) && fn < f
            w = wn;  F = Fn;  f = fn;  lam = max(lam/5,1e-15);  ok = true;  break %#ok<NASGU> % F reused unless the outer loop exits
        end
        lam = lam*8;
    end
    if ~ok || f < tol, break; end
    [F,J] = bm_resid(w,P,rv);
end
end
