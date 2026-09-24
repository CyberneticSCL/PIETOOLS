function X = pielr_evalop(prog,X)                                           % CC, 09/23/2026
% PIELR_EVALOP  Substitute a decision vector into ANY LPI decision object,
% dispatching on its class: dopvar (1-D, 6 leaves), dopvar2d (2-D, 36 leaves),
% or a bare dpvar (a scalar decision variable such as the l2gain gamma).
%
% WHY.  Verification has to push the candidate back through the real operators
% in BOTH dimensions and for programs that carry more than a Gram: the l2gain
% LPI substitutes Rop (Gram), Qop (a free lpivar operator) and gam (a scalar)
% before its residual can be formed.  opcheck_2d could reach only dopvar2d via
% pielr_evalgram, and the 1-D checker used lpigetsol, whose dpvar2poly path is
% the one measured not to finish in 35 minutes on a ~277k-dvar operator
% (pielr_evalgram's header).  One dispatcher keeps both dimensions on the fast
% sparse path and keeps the two from drifting.
%
% An already-substituted object (opvar / opvar2d / double) passes through.

if isa(X,'opvar') || isa(X,'opvar2d') || isa(X,'double') || isa(X,'polynomial')
    return                                  % nothing to substitute
end
dtab = prog.decvartable;
RRx  = prog.solinfo.RRx;

if isa(X,'dopvar2d')
    X = pielr_evalgram(prog,X);   return
end
if isa(X,'dpvar')
    % scalar / matrix decision variable with no operator structure
    X = pielr_evaldp(X,dtab,RRx);  return
end
if ~isa(X,'dopvar')
    error('pielr_evalop:input', ...
          'cannot substitute into an object of class %s',class(X));
end
% ---- 1-D dopvar: the 6 leaves P, Q1, Q2, R.R0, R.R1, R.R2 ----------------
% Enumeration taken from the dopvar classdef properties (P, Q1, Q2, R struct
% of three); pielr_maxop asserts the count on the result, so a class that
% grows a parameter fails loudly instead of being silently skipped.
%
% A fresh opvar is populated cell by cell rather than opvar(X)'d: unlike
% opvar2d, the 1-D constructor does not accept a dopvar (it reads a char
% argument list and errors with "Input must be strings").  This mirrors
% getsol_lpivar, which builds its result the same way and for the same reason.
Y = opvar();
Y.I = X.I;   Y.var1 = X.var1;   Y.var2 = X.var2;   Y.dim = X.dim;
Y.P  = pielr_evaldp(X.P ,dtab,RRx);
Y.Q1 = pielr_evaldp(X.Q1,dtab,RRx);
Y.Q2 = pielr_evaldp(X.Q2,dtab,RRx);
Y.R.R0 = pielr_evaldp(X.R.R0,dtab,RRx);
Y.R.R1 = pielr_evaldp(X.R.R1,dtab,RRx);
Y.R.R2 = pielr_evaldp(X.R.R2,dtab,RRx);
Y.P = double(Y.P);                % getsol_lpivar's own last line
X = Y;
end
