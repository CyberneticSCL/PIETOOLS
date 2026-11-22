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
%
% CC, 11/22/2026: the dopvar and dopvar2d branches deleted (unreachable since
%   the opvar classes left the path, f83486c8; pielr_evalgram retired with
%   pielr_certify, maintainer: "port or retire"). Containers and dpvar remain.

% BEGIN CC, 11/08/2026: container operators.  A 'copvar' or 'sopvar' holds
% no decision variable; a 'cdopvar' or 'sdopvar' is substituted by
% GETSOL_LPIVAR_SOP, which reads prog.solinfo.RRx against prog.decvartable --
% the vector pielr_opcheck writes there before calling this.  Checked first:
% the opvar classes below left the path (f83486c8) and their isa tests now
% return false rather than error.
if isa(X,'copvar') || isa(X,'sopvar')
    return                                  % nothing to substitute
end
if isa(X,'cdopvar') || isa(X,'sdopvar')
    X = getsol_lpivar_sop(prog,X);   return
end
% END CC, 11/08/2026
if isa(X,'opvar') || isa(X,'opvar2d') || isa(X,'double') || isa(X,'polynomial')
    return                                  % nothing to substitute
end
dtab = prog.decvartable;
RRx  = prog.solinfo.RRx;

% BEGIN CC, 11/22/2026: deleted, both CC 09/23/2026: the dopvar2d branch
% (pielr_evalgram, retired to orphaned_file/private/) and the 1-D dopvar
% branch (an opvar() filled leaf by leaf through pielr_evaldp). Their
% classes left the path (f83486c8): isa is false, neither was reachable.
if isa(X,'dpvar')
    % scalar / matrix decision variable with no operator structure
    X = pielr_evaldp(X,dtab,RRx);  return
end
error('pielr_evalop:input', ...
      'cannot substitute into an object of class %s',class(X));
% END CC, 11/22/2026
end
