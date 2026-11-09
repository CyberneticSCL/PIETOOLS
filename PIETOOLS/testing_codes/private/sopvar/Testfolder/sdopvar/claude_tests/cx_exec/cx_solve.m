function r = cx_solve(prog,sopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = CX_SOLVE(PROG,SOPTS) solves an assembled LPI program with lpisolve and
% returns a CERTIFIED verdict plus the SDP shape and timings, identically for
% the container and the stock programs.
%
% OUTPUT struct R:
%   st      +1 certified feasible (pinf = 0, numerr = 0), -1 certified
%           infeasible (pinf = 1, numerr = 0), 0 otherwise: a verdict with a
%           numerical error certifies nothing, so is never counted as either;
%   pinf, dinf, numerr, feasratio   sossolve's info fields (for Mosek,
%           feasratio is MSK_DINF_INTPNT_OPT_STATUS, sossolve.m:368);
%   cpusec  the solver's own clock (sossolve.m:366 for Mosek);
%   wall    wall time of lpisolve, including sossolve's assembly;
%   obj     the objective value when the program has one (prog.objective'*x);
%   shape   cx_shape of the solved program (post-solve, so the slack blocks
%           sossolve adds for 'ineq' constraints are included);
%   sol     the solved program.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = tic;
evalc('sol = lpisolve(prog,sopts);');
wall = toc(t);
q = struct();
if isfield(sol,'solinfo') && isfield(sol.solinfo,'info'),   q = sol.solinfo.info;  end
g = @(f) getf(q,f);
r = struct();
r.pinf = g('pinf');     r.dinf = g('dinf');     r.numerr = g('numerr');
r.feasratio = g('feasratio');   r.cpusec = g('cpusec');     r.wall = wall;
r.st = 0;
if r.numerr==0 && r.pinf==0,    r.st = +1;  end
if r.numerr==0 && r.pinf==1,    r.st = -1;  end
r.obj = NaN;
if isfield(sol,'objective') && any(sol.objective) && isfield(sol.solinfo,'RRx') ...
        && ~isempty(sol.solinfo.RRx)
    nob = numel(sol.objective);
    r.obj = full(sol.objective(:)'*sol.solinfo.RRx(1:nob));
end
r.shape = cx_shape(sol);
[r.rel_b,r.psd_min,r.psd_relmin,r.trivial] = cx_resid(sol);
r.sol = sol;
end


function v = getf(q,f)
if isfield(q,f) && ~isempty(q.(f)),     v = double(q.(f));  else,   v = NaN;    end
end
