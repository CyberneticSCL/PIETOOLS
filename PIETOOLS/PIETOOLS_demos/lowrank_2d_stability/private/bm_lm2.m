function [w,f,it,why] = bm_lm2(P,rv,w0,maxit,tol,rtol,win)                  % CC, 09/22/2026
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
%
% CC, 09/22/2026: add a WINDOWED STAGNATION TEST, and return why we stopped.
%   Before this, the only exits were a stalled damping loop and f < tol with
%   tol = 1e-16 from the caller, which never fires -- so maxit was a HARD CUT,
%   not convergence.  MEASURED (1-D heavy, m=56, rank 2/block, chunked runs to
%   4000 iterations at lam/lam* = 0.20/0.40/0.60/0.80/0.95): the residual is
%   STILL DESCENDING at 4000 at every lambda, improving 25x to 117x past
%   iteration 400, with log-slope -0.07 to -0.30 log10 per 1000 iterations and
%   NO plateau anywhere -- the slope STEEPENS with lambda.  Every reach/limit
%   this package has reported was therefore measuring where the cut landed.
%   The fix is a convergence criterion, not a bigger constant: with it, a large
%   maxit becomes safe, because easy points exit early (they were burning the
%   full budget for a residual already 25x below the gate) while hard points
%   run until they genuinely flatten.
%   THRESHOLD CHOICE, from those slopes: the flattest measured descent is
%   -0.07 log10/1000 its = 1.6% per 100 its, so a 0.1%-per-100-its floor sits
%   ~16x below the flattest real progress and cannot stop a descending run.
%   MEASURED CONSEQUENCE, worth stating plainly: IN 1-D the stagnation test
%   NEVER FIRES -- re-entering at an apparently converged point still improves
%   (7.58e-08 -> 7.21e-08 over another 4000 its) at all five lambdas.  There it
%   is a safety net for a genuinely flat problem, not the exit that matters.
%   CC, 09/22/2026: THE SENTENCE ABOVE WAS WRITTEN "on this problem family" AND
%   THAT WAS TOO BROAD -- it was measured in 1-D only, and IN 2-D THE TEST DOES
%   FIRE.  On the anisotropic 2-D heat variant at rank [2 2], seeds 11 and 22
%   return BIT-IDENTICAL raw and op values at lmit 2000 and 8000 (3.352e-06 /
%   1.166e-06 and 3.301e-06 / 1.307e-06), which is only possible if the run
%   exited early at the same point both times; 'tol' cannot be it, since raw
%   3.3e-06 never crosses lmtol 1e-6.  Seed 33 on the same problem keeps
%   descending (4.888e-06 -> 1.082e-06 -> 1.000e-06), so the test is also not
%   firing spuriously on a run that still has progress left.  This is what
%   makes a LARGE maxit safe rather than wasteful, and it is the argument for
%   raising the default: in 2-D the budget is declined when it is not needed.  THE EXIT THAT MATTERS IS 'tol': the caller does not
%   need ||F|| -> 0, it needs the OPERATOR GATE, and measured on these programs
%   op_rel ~ 0.35-0.40 * raw_rel (raw 1.948e-06 -> op 7.606e-07; raw 2.639e-06
%   -> op 9.291e-07).  So a raw target at the gate value leaves ~2.5x margin.
%   Passing tol = 1e-16, as this was called with, disables that exit entirely
%   and is what turned maxit into the operative stopping rule.
w = w0;  lam = 1e-8;
if nargin<4||isempty(maxit), maxit=200; end
if nargin<5||isempty(tol),   tol=1e-13;  end
if nargin<6||isempty(rtol),  rtol=1e-3;  end  % min relative decrease/window  % CC, 09/22/2026
if nargin<7||isempty(win),   win=100;    end  % window, iterations            % CC, 09/22/2026
why = 'maxit';                                % overwritten by every other exit % CC, 09/22/2026
[F,J] = bm_resid(w,P,rv);  f = norm(F);
fwin = f;                                     % residual at window start      % CC, 09/22/2026
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
    % CC, 09/22/2026: DAMPING POLICY LEFT AS SHIPPED, on measurement.  A
    % gain-ratio (Nielsen) update was implemented here and REVERTED: it is
    % WORSE.  Measured, 1-D heavy, rank 2/block, identical start, 300 its,
    % ||F|| at the end relative to this policy (>1 = better):
    %   policy (accept factor, reject factor)   lam/lam*=0.20   0.95
    %   /5, x8 flat            [SHIPPED]            1.00         1.00
    %   /5, x2 escalating                           0.83         0.93
    %   /10, x2 escalating                          0.52         0.87
    %   Nielsen max(1/3,1-(2rho-1)^3), x2 esc       0.73         0.87
    %   /5, x2 flat                                 0.73         0.85
    % WHY the intuition was wrong: at the measured gain ratios (median rho
    % 0.65-0.76) the Nielsen factor is 0.85-0.97, i.e. NEARLY NEUTRAL -- it
    % only crashes lambda as rho -> 1.  The flat /5 = 0.2 is already far more
    % aggressive there.  And a gentler REJECTION factor backfires: x2 does not
    % escape a bad region in one step, so rejections rise from 64% to 71-95%
    % of iterations and tries/iteration from 2 to 3.  The single x8 jump is
    % doing real work.  Do not re-tune this without re-running that sweep
    % (scratchpad/reach1d/q6.m) -- the shipped values won five ways.
%   if ~ok || f < tol, break; end                                           % CC, 09/22/2026 (was)
    if ~ok,      why = 'damping';  break; end   % 40 tries found no decrease % CC, 09/22/2026
    if f < tol,  why = 'tol';      break; end                               % CC, 09/22/2026
    % Windowed stagnation: only checked every 'win' iterations, because a
    % per-iteration test would fire on the ordinary sawtooth of a damped
    % Gauss-Newton step and cut a descending run short -- the defect this is
    % here to prevent, not to reproduce.
    if mod(it,win)==0                                                       % CC, 09/22/2026
        if (fwin - f) <= rtol*fwin, why = 'stagnant';  break; end           % CC, 09/22/2026
        fwin = f;                                                           % CC, 09/22/2026
    end                                                                     % CC, 09/22/2026
%   [F,J] = bm_resid(w,P,rv);                                               % CC, 09/22/2026 (was)
    % F is already Fn from the accepted trial; ask for J only.  See bm_resid's
    % header for the measured split that makes this ~1/6 of a step.
    [~,J] = bm_resid(w,P,rv,true);                                          % CC, 09/22/2026
end
end
