% DEMO_YOUR_SYSTEM  Worked template: certify YOUR OWN 2-D system.
%
% This is a SCRIPT: everything it computes -- cert, cert2, cert3, PIE, ... --
% stays in your workspace afterwards, so you can keep testing on the results.
%
% Two ways in, both shown below and both live:
%   (A) define a PDE with pde_var, convert() it, hand the PIE to pielr_certify;
%   (B) hand pielr_certify RAW opvar2d OPERATORS T and A directly -- no
%       convert() call -- by assembling a bare struct.
%
% The example PDE is deliberately NOT the demo's heat equation: an ANISOTROPIC
% diffusion with a reaction term,
%     x_t = x_{s1s1} + 0.5*x_{s2s2} + lam*x   on (0,1)^2, Dirichlet edges.
% Separation of variables gives the spectrum lam - (j^2 + 0.5*k^2)*pi^2,
% j,k >= 1, so it is exponentially stable iff lam < 1.5*pi^2 = 14.804.
% lam = 2 below is well inside.
%
% MEASURED on this exact example (fresh run): the cold call certifies at
% r = [3 3], op rel 9.9e-07, discovery ~25 min.  Rank 2 was NOT reached --
% three BM seeds landed at 1.4e-06..4.9e-06, just above the gate -- so [3 3]
% is an honest upper bound, not a floor.  Expect discovery time and the
% achieved rank to vary with the system; that variation is the point.
%
% CC, 09/19/2026: converted from a function to a SCRIPT (maintainer request)
%                 so the certificates persist in the workspace; the
%                 maintainer's follow-on experiment (a nearby lam, warm-started
%                 from the previous certificate) kept as section (5).

here = fileparts(mfilename('fullpath'));
cd(here);              % 2-D conversion measured to fail from a cluttered cwd
pielr_path();

% ============ (1) define YOUR PDE with pde_var ===============================
% Replace this block with your own dynamics / boundary conditions.
lam = 2;
pvar s1 s2
clear stateNameGenerator        % pde_var keeps counters between definitions
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
sys = [diff(x,'t')==diff(x,s1,2)+0.5*diff(x,s2,2)+lam*x;
       subs(x,s1,0)==0;  subs(x,s1,1)==0;
       subs(x,s2,0)==0;  subs(x,s2,1)==0];

% ============ (2) convert to a PIE ===========================================
PIE = convert(sys,'pie');

% ============ (3) certify, cold ==============================================
% No face is supplied, so pielr_certify searches for one (discovery) and then
% certifies it.  Expect minutes, dominated by discovery -- the report prints
% setup / discovery / certification separately.  All options are defaulted;
% see 'help pielr_certify' for the knobs (rank to try first, search ceiling,
% residual gate, settings/degrees, BM vs min-trace route).
cert = pielr_certify(PIE);
assert(cert.ok,'demo_your_system: the cold call did not certify -- see the report above.');

% ============ (4) OR: hand in your own T and A operators directly ============
% If you already have the PIE operators as opvar2d objects -- from your own
% construction, not from convert() -- assemble a bare struct.  Only T and A are
% used (autonomous stability); vars and dom are optional and default from T
% itself, so the two lines marked (optional) can be dropped.
Top = PIE.T;                    % <-- your own opvar2d goes here
Aop = PIE.A;                    % <-- your own opvar2d goes here
raw = struct();
raw.T    = Top;
raw.A    = Aop;
raw.vars = [Top.var1, Top.var2];   % (optional) primary/dummy variables, 2x2
raw.dom  = Top.I;                  % (optional) spatial domain, 2x2
% The certificate found in step (3) is passed straight back in, so this second
% call is certification only (seconds, no discovery).  For a genuinely
% different system, drop the second argument and let it search.
cert2 = pielr_certify(raw,cert);       % a previous cert IS a valid input
assert(cert2.ok,'demo_your_system: the raw-operators call did not certify.');

fprintf(['\nBoth entry points certified: rank %s (converted PIE, cold) and ' ...
         'rank %s (raw T/A struct).\n'],mat2str(cert.r),mat2str(cert2.r));

% ============ (5) why faces do NOT transfer across physics ===================
% Demonstration, not a workflow: change the system (here lam 2 -> 1.9) and
% try re-certifying from the previous certificate.  THIS FAILS, and that is
% the point -- a face is specific to the operating point (MEASURED: the
% lam=2 face gives op rel = 1 at lam=1.9 here, op rel 1.7e-2 on the heat
% demo, while a cold call re-certifies lam=1.9 at the SAME rank).  Reuse a
% certificate only on the SAME physics (saved certs, replicated states); a
% changed system pays its own discovery.  The failure is honest -- the
% restriction can only lose feasibility, never fake a certificate -- and
% says NOTHING about the new system's stability.
% CC, 09/19/2026: assert removed -- the failure is an expected, informative
%                 outcome, not a script error (measured above).
% CC, 09/20/2026: reframed as a demonstration of non-transfer; cross-physics
%                 face reuse is abandoned as a workflow (maintainer decision).
pvar s1 s2
clear stateNameGenerator
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
sys = [diff(x,'t')==diff(x,s1,2)+0.5*diff(x,s2,2)+1.9*x;
       subs(x,s1,0)==0;  subs(x,s1,1)==0;
       subs(x,s2,0)==0;  subs(x,s2,1)==0];
PIE19 = convert(sys,'pie');
cert3 = pielr_certify(PIE19,cert2);
if cert3.ok
    fprintf('Nearby system (lam = 1.9) re-certified from the old face in %.1f s.\n', ...
            cert3.t_certify);
else
    fprintf(['Nearby system (lam = 1.9): the old face contains NO certificate ' ...
             '(op rel %.3g)\n-- the expected face-transfer failure, and not a ' ...
             'statement about stability.\nA cold call (cert3 = ' ...
             'pielr_certify(PIE19);) re-certifies at the same rank in minutes.\n'], ...
            cert3.op_rel);
end
fprintf('DEMO_YOUR_SYSTEM_DONE\n');
