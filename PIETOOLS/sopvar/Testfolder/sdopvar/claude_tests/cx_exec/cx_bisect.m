function B = cx_bisect(f,g0,rtol,maxit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B = CX_BISECT(F,G0,RTOL,MAXIT) locates the feasibility threshold of a
% fixed-gamma LPI on CERTIFIED verdicts only. F(g) returns a cx_solve
% struct (field st: +1 feasible, -1 infeasible, 0 uncertified).
%
% From G0 it doubles up to a certified-feasible gamma and halves down to a
% certified-infeasible one, then bisects geometrically until hi/lo - 1 <=
% RTOL (default 1e-3). An uncertified midpoint is recorded and bisection
% stops there: [lo,hi] is then the certified bracket, and B.gap says the
% threshold was not resolved further. Monotone in gamma by construction of
% the Hinf/H2 LPIs (a larger gamma relaxes the constraint).
%
% OUTPUT struct: lo, hi (certified bracket; NaN if not found), est = sqrt(lo*hi),
% resolved (true if hi/lo-1 <= rtol), trace (rows [g st]), nsolve, t (wall s).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || isempty(rtol),   rtol = 1e-3;    end
if nargin<4 || isempty(maxit),  maxit = 40;     end
tr = zeros(0,2);    t0 = tic;
ev = @(g) f(g).st;
lo = NaN;   hi = NaN;
g = g0;     s = ev(g);  tr(end+1,:) = [g s];
if s==+1,       hi = g;
elseif s==-1,   lo = g;
end
% Expand to a certified bracket.
k = 0;
while isnan(hi) && k<12
    g = 2*max([g, lo]);     s = ev(g);  tr(end+1,:) = [g s];     k = k+1;
    if s==+1,   hi = g;     elseif s==-1,   lo = g;     end
end
k = 0;
while isnan(lo) && ~isnan(hi) && k<12
    g = min([g, hi])/2;     s = ev(g);  tr(end+1,:) = [g s];     k = k+1;
    if s==-1,   lo = g;     elseif s==+1,   hi = g;     end
end
% Bisect.
it = 0;     stopped = false;
while ~isnan(lo) && ~isnan(hi) && hi/lo-1 > rtol && it<maxit
    g = sqrt(lo*hi);    s = ev(g);  tr(end+1,:) = [g s];     it = it+1;
    if s==+1,       hi = g;
    elseif s==-1,   lo = g;
    else,           stopped = true;     break
    end
end
B = struct('lo',lo,'hi',hi,'est',sqrt(lo*hi),'resolved',~stopped && hi/lo-1<=rtol, ...
           'trace',tr,'nsolve',size(tr,1),'t',toc(t0));
end
