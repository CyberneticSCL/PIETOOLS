% ab_bisect: lb_bisect and its helpers of lb_run.m (11/17 harness), extracted
% for the A/B scripts; the first function renamed. MMP, 11/17/2026.
function out = ab_bisect(c, S0, opts)
% Largest certified parameter value in [c.lo, c.hi]: lo must certify, hi
% must not; then halve until hi-lo <= tol or opts.maxsolves solves.
lo = c.lo;  hi = c.hi;
if isfinite(c.docnum) && c.docnum~=0
    tol = opts.reltol*abs(c.docnum);
else
    tol = opts.reltol*(hi-lo);
end
steps = {};  n = 0;  ndec = NaN;  last = [];
[ok, r, e] = lb_step(c, S0, lo);  n = n+1;  steps{end+1} = lb_stp(lo,r,e);
if ~isempty(r), ndec = r.ndec; last = r; end
if ~isempty(e) && isempty(r)
    rethrow(e);                      % the default bracket itself fails to run
end
if ~ok
    out = lb_bout('lo_fail', NaN, n, ndec, steps, last, c.lo, c.hi, lo, hi);  return
end
[ok, r, e] = lb_step(c, S0, hi);  n = n+1;  steps{end+1} = lb_stp(hi,r,e);
if ~isempty(r), ndec = r.ndec; last = r; end
if ok
    out = lb_bout('hi_pass', hi, n, ndec, steps, last, c.lo, c.hi, lo, hi);  return
end
while hi-lo > tol && n < opts.maxsolves
    mid = (lo+hi)/2;
    [ok, r, e] = lb_step(c, S0, mid);  n = n+1;  steps{end+1} = lb_stp(mid,r,e); %#ok<AGROW>
    if ~isempty(r), ndec = r.ndec; last = r; end
    if ok, lo = mid; else, hi = mid; end
end
if hi-lo <= tol, st = 'converged'; else, st = 'max_solves'; end
out = lb_bout(st, lo, n, ndec, steps, last, c.lo, c.hi, lo, hi);
end

function [ok, r, e] = lb_step(c, S0, v)
% One bisection solve at parameter value v; an error counts as a failure.
r = [];  e = [];
try
    PIE = lb_build_system(c, v, true);
    r = lb_exec(c.exec, PIE, S0);
    ok = strcmp(r.status,'certified');
catch e
    ok = false;
end
end

function s = lb_stp(v, r, e)
code = struct('certified','C','inaccurate','I','infeasible','F','unknown','U','nosolinfo','N');
if ~isempty(e) || isempty(r), s = sprintf('%.6g:E',v);
else, s = sprintf('%.6g:%s(%.0fs)',v,code.(r.status),r.sec); end
end

function out = lb_bout(st, res, n, ndec, steps, last, lo0, hi0, lo, hi)
out.status = st;  out.result = res;  out.n = n;  out.ndec = ndec;
out.note = sprintf('bracket [%.6g %.6g] final [%.6g %.6g] | steps %s', lo0, hi0, lo, hi, strjoin(steps,' '));
if ~isempty(last), out.note = [out.note,' | last: ',lb_rnote(last)]; end
end

function s = lb_rnote(r)
s = sprintf('fr=%.3g pinf=%g dinf=%g numerr=%g iter=%g cpu=%.1fs m=%g nx=%g', ...
    r.fr, r.pinf, r.dinf, r.numerr, r.iter, r.cpusec, r.m, r.nx);
if isfield(r,'called') && startsWith(r.called,'PIETOOLS_') && contains(r.called,'_sop')
    s = [r.called,' ',s];
end
end
