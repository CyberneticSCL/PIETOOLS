function T = cx_stability_thresh(c,rtol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = CX_STABILITY_THRESH(C,RTOL) bisects the destabilizing plant parameter
% frac = C.plant{2} of a stability case (cx_cases_stability) on CERTIFIED
% verdicts (cx_solve st), on both paths: the stock executive's own program
% (cx_stock_capture, then solved by cx_solve) and the container
% transcription cx_<exec>.
%
% Why not cx_compare's thresh: it calls cx_bisect, which assumes the LPI is
% feasible at LARGE parameter values (true of gamma). For 'rd' a larger frac
% is less stable (u_t = u_ss + frac*pi^2 u, threshold frac = 1), so that
% bisection halves toward 0 looking for an infeasible point it never finds.
% Here feasibility is at SMALL frac: from frac0 (must be certified feasible)
% frac grows by x1.5 until certified infeasible, then geometric bisection
% until hi/lo-1 <= RTOL (default 1e-2); the first uncertified midpoint stops
% it and [lo,hi] stays the certified bracket (resolved = false).
%
% OUTPUT struct: stock, cx, each with lo (feasible), hi (infeasible),
% resolved, trace (rows [frac st]), nsolve, t; agree (brackets overlap).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(rtol),   rtol = 1e-2;    end
if ~isfield(c,'solver') || isempty(c.solver),   c.solver = 'mosek';   end
st = cx_settings(c.setname,c.solver);   sopts = st.sos_opts;
pl = c.plant;   mkP = @(fr) cx_plant(pl{1},fr,pl{3:end});
fs = @(fr) cx_solve(cx_stock_capture(c.exec,mkP(fr),st),sopts);
fc = @(fr) cx_solve(feval(['cx_' c.exec],mkP(fr),st),sopts);
T.stock = bis(fs,pl{2},rtol);
T.cx = bis(fc,pl{2},rtol);
A = T.stock;    B = T.cx;
T.agree = all(isfinite([A.lo A.hi B.lo B.hi])) && A.lo < B.hi && B.lo < A.hi;
end

function B = bis(f,f0,rtol)
t0 = tic;   tr = zeros(0,2);
lo = NaN;   hi = NaN;   stopped = false;
fr = f0;    s = f(fr).st;   tr(end+1,:) = [fr s];
if s==+1,   lo = fr;    elseif s==-1,   hi = fr;    end
k = 0;
while isnan(hi) && ~isnan(lo) && k<12              % grow to an infeasible point
    fr = 1.5*fr;    s = f(fr).st;   tr(end+1,:) = [fr s];   k = k+1;
    if s==+1,   lo = fr;    elseif s==-1,   hi = fr;    end
end
k = 0;
while isnan(lo) && ~isnan(hi) && k<12              % shrink to a feasible point
    fr = fr/1.5;    s = f(fr).st;   tr(end+1,:) = [fr s];   k = k+1;
    if s==+1,   lo = fr;    elseif s==-1,   hi = fr;    end
end
while ~isnan(lo) && ~isnan(hi) && hi/lo-1 > rtol
    fr = sqrt(lo*hi);   s = f(fr).st;   tr(end+1,:) = [fr s];
    if s==+1,       lo = fr;
    elseif s==-1,   hi = fr;
    else,           stopped = true;     break
    end
end
B = struct('lo',lo,'hi',hi,'resolved',~stopped && hi/lo-1<=rtol,'trace',tr, ...
           'nsolve',size(tr,1),'t',toc(t0));
end
