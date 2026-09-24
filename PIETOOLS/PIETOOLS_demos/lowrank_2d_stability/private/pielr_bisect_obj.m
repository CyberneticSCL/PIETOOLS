function [V,R,q,why,notes,P] = pielr_bisect_obj(prog,H,D,A,opts,vb)        % CC, 09/23/2026
% PIELR_BISECT_OBJ  Certify an LPI that carries a scalar objective, by
% bisection on the objective variable.
%
% WHY BISECTION AND NOT "JUST SOLVE IT".  Nothing in the low-rank pipeline
% optimises anything.  bm_lm2 minimises ||W(A(X)-b)||, in which the objective
% vector c never appears, and restrict_solve pins s by a determined solve --
% np == rM at every rank measured, so on a face there is not one degree of
% freedom left to spend on the objective.  Run as-is on an l2gain LPI the
% pipeline therefore returns whatever gamma the face happens to carry:
% MEASURED on Ex_Transport_Eq_with_Disturbance at 'light', gamma = 1.42095
% from the low-rank search against 0.516333 from the interior-point solve of
% the same program -- 2.75x, sound as an upper bound and useless as a gain.
%
% Fixing gamma turns the optimisation LPI back into the feasibility LPI the
% machinery is built for, and bisecting recovers the optimum.  Each certified
% gamma is a genuine UPPER bound, because a face can only lose feasible
% points, so every iterate of the bisection is a valid answer and the
% bisection only tightens it.
%
% HOW gamma IS FIXED.  The objective coordinate jc is pinned by appending one
% equality row e_jc' x = g to the SeDuMi data.  That is a column of At, since
% the constraint is At'*x = b.  Pinning rather than substituting keeps the
% coordinate in place, so the layout, the operator handles and the gate are
% all untouched.
%
% COST, which is the real objection to this design.  Each trial rebuilds
% bm_setup, and with pre=1 that is a DENSE m x m eigendecomposition -- O(m^3)
% time and O(m^2) memory per bisection step.  At 1-D sizes (m ~ 110) it is
% free; at the 2-D sizes this package targets (m ~ 6000) it is the dominant
% cost and a bisection of 12 steps multiplies it by 12.  The whitener is
% computed from Ssym*Ssym' and only ONE COLUMN of At changes between trials,
% so this is a rank-one update in principle; left unexploited here and noted
% because it is the first thing to fix if the 2-D arm is wanted.
%
% OUTPUT  V,R,q,why,notes as pielr_discover;  P the bm_setup package of the
%         ACCEPTED trial, which the caller needs to lift the face.

notes = {};  V = [];  R = [];  q = [];  why = [];  P = [];
jc = find(D.c);
if numel(jc)~=1
    error('pielr_bisect_obj:obj', ...
        ['bisection handles a SCALAR objective coordinate; this program has ' ...
         '%d nonzero entries in c.  A vector objective needs a different ' ...
         'scheme (scalarise, or put the objective into the restricted solve).'], ...
         numel(jc));
end
sgn = sign(D.c(jc));            % +1: minimising, which is every case here
if sgn <= 0
    error('pielr_bisect_obj:sign','objective coefficient is not positive; expected a minimisation');
end

lo = 0;   hi = [];   bestall = struct('V',[],'R',[],'q',[],'P',[],'g',inf);
allwhy = [];
nit = 0;

% ---- 1. unpinned run, to bracket ---------------------------------------
if vb, fprintf('  [gamma] unpinned run, to bracket ...\n'); end
[V0,R0,q0,w0why,n0,P0] = trial(prog,H,D,A,opts,vb,[],[]);
allwhy = [allwhy w0why];  notes = [notes n0];
if ~isempty(R0) && R0.ok && isfield(R0,'aux') && isfield(R0.aux,'gam')
    hi = R0.aux.gam;
    bestall = struct('V',{V0},'R',R0,'q',q0,'P',P0,'g',hi);
    if vb, fprintf('  [gamma] bracket from the unpinned run: %.6g\n',hi); end
else
    % no bracket yet: climb a geometric ladder until something certifies
    g = 1;
    for k = 1:12
        nit = nit+1;
        [Vk,Rk,qk,wk,nk,Pk] = trial(prog,H,D,A,opts,vb,jc,g);
        allwhy = [allwhy wk];  notes = [notes nk];   %#ok<AGROW>
        if vb, fprintf('  [gamma] ladder g = %-10.6g -> %s\n',g,verd(Rk)); end
        if ~isempty(Rk) && Rk.ok
            hi = g;  bestall = struct('V',{Vk},'R',Rk,'q',qk,'P',Pk,'g',g);
            break
        end
        lo = g;  g = g*2;
    end
end
if isempty(hi)
    notes{end+1} = ['no gamma certified on the bracketing ladder up to 2^12; ' ...
        'reported as NOT REACHED, never as a proven bound'];
    why = allwhy;   return
end

% ---- 2. bisect ----------------------------------------------------------
% relative tolerance on gamma; 1e-3 puts the answer well inside the spread
% between settings levels, which is what the reference values differ by.
rtol = 1e-3;   maxit = 20;
while (hi-lo) > rtol*max(hi,eps) && nit < maxit
    nit = nit+1;
    g = 0.5*(lo+hi);
    [Vk,Rk,qk,wk,nk,Pk] = trial(prog,H,D,A,opts,vb,jc,g);
    allwhy = [allwhy wk];  notes = [notes nk];       %#ok<AGROW>
    if vb, fprintf('  [gamma] bisect g = %-10.6g -> %s\n',g,verd(Rk)); end
    if ~isempty(Rk) && Rk.ok
        hi = g;  bestall = struct('V',{Vk},'R',Rk,'q',qk,'P',Pk,'g',g);
    else
        lo = g;
    end
end
V = bestall.V;  R = bestall.R;  q = bestall.q;  P = bestall.P;
why = allwhy;
notes{end+1} = sprintf( ...
    ['gamma %.6g by bisection over %d trials (bracket [%.4g, %.4g], rel tol %g); ' ...
     'each certified gamma is an UPPER bound, so this is one too'], ...
    bestall.g,nit,lo,hi,rtol);
end

% =========================================================================
function [V,R,q,why,notes,P] = trial(prog,H,D,A,opts,vb,jc,g)
% One feasibility run with the objective coordinate pinned (or not, when jc
% is empty).  Verbosity is dropped to the per-attempt table only at the top
% level: a bisection prints 15 tables otherwise.
Dt = D;
if ~isempty(jc)
    e = sparse(jc,1,1,D.L.Ntot,1);
    Dt.At = [D.At, e];
    Dt.b  = [D.b;  g];
end
P = bm_setup(Dt.At,Dt.b,Dt.L.N,Dt.L.Kf,1,Dt.L);
[V,R,q,why,notes] = pielr_discover(prog,H,Dt,P,A,opts,vb && isempty(jc));
end

function s = verd(R)
if isempty(R), s = 'no certificate';
elseif R.ok,   s = sprintf('CERTIFIES (rel %.3e)',R.rel);
else,          s = sprintf('fail (rel %.3e)',R.rel);
end
end
