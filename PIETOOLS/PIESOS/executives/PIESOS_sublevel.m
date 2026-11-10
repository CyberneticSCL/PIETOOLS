function [prog_sol,info] = PIESOS_sublevel(PIE,V,g,h,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG_SOL,INFO] = PIESOS_sublevel(PIE,V,G,H,SETTINGS) finds the largest
% sublevel set {V <= c} of the functional V that lies in the semialgebraic
% set
%   {g_i(x) >= 0, i=1,...,ng,  h_j(x) = 0, j=1,...,nh}
% of the PIE state x=PIE.x, in the sense that
%   V(x) <= c, h_j(x) = 0 for all j     implies     g_i(x) >= 0 for all i.
% The equality constraints h_j describe the states considered (for example,
% the ball N(u) <= R^2 of PIESOS_stability_local is g = R^2 - N(u), without
% h). If V is a Lyapunov functional for which the derivative dV along the
% PIE is nonpositive on {g >= 0, h = 0}, for example the output V of
% PIESOS_stability_local with g = R^2 - N(u), then {V <= c} is an invariant
% set in this set, and thus a region of attraction.
% With constants a_i = g_i(0) >= 0, the inclusion for g_i is verified using
% a multiplier p_i >= 0 on V - c, and free multipliers q_ij on h_j, by
%   a_i - p_i*c >= 0,    (g_i - a_i) + p_i*V + sum_j q_ij*h_j >= 0 (SOS).
% Then g_i + p_i*(V - c) + sum_j q_ij*h_j >= 0, so that g_i >= 0 on
% {V <= c, h = 0}. The smallest p_i for which the second condition holds is
% found by one SOS program per inequality, giving c_i = a_i/p_i, and the
% largest sublevel set has c = min_i c_i.
%
% INPUT
% - PIE:        'npie_struct' object representing the PIE, as output by
%               `convert` for a nonlinear PDE. A 'pde_struct' object is
%               converted to a PIE;
% - V:          'polyopvar' object in the PIE state x=PIE.x with V(0)=0,
%               e.g. info.V as output by PIESOS_stability_local;
% - g:          scalar 'polyopvar' object in x, or cell of such objects,
%               defining the inequalities g_i >= 0, with g_i(0) >= 0. The
%               constant term is allowed, e.g. g = R^2 - weighted_sobolev_
%               norm(PIE.T*PIE.x,alpha);
% - h:          (optional) scalar 'polyopvar' object in x, or cell of such
%               objects, defining the equalities h_j = 0, with h_j(0)=0.
%               Defaults to {};
% - settings:   (optional) 'struct' of settings, as output by
%               `PIESOS_settings`. Only the following fields are used:
%                 LF.mon_deg,           degree of the kernels in V, to be set
%                                       to the value used to obtain V;
%                 sos.mon_deg(3),       degree of the kernels in the SOS
%                 sos.psatz(3)          constraint, the larger of
%                                       sos.mon_deg(3) and LF.mon_deg+2 (a
%                                       lower degree, or one chosen
%                                       automatically, was infeasible for
%                                       the V of PIESOS_stability_local), and
%                                       psatz multipliers;
%                 solver, simplify, tol_feas, verbose
%               Defaults to PIESOS_settings();
%
% OUTPUT
% - prog_sol:   PIESOS program structure of the last solved program;
% - info:       'struct' with fields
%       feasible:   logical, true if all programs were solved, so that the
%                   sublevel set {V <= c} is verified;
%       c:          scalar, the largest value verified, min_i c_i (0 if not
%                   all programs were solved). Inf if g_i >= 0 holds
%                   everywhere on {h = 0} for all i;
%       c_i:        ng x 1 array, the values a_i/p_i (NaN if not solved);
%       time:       the time taken, in seconds;
%
% NOTES
% - The multipliers on V - c are constants, so that c is exact if V and g_i
%   are quadratic (S-lemma), but otherwise conservative: only the quadratic
%   part of g_i - a_i is dominated by V. A multiplier that is a polynomial
%   in the states would require SOS conditions on products such as ||u||^4,
%   which piesos_ineq did not certify for PDE states (even at high degrees).
% - The region {V <= c} is only a region of attraction if dV <= 0 on
%   {g >= 0, h = 0}. This is not checked here.
% - For radius of the region of attraction in a norm N with V <= C*N on the
%   ball, e.g. C=info.C of PIESOS_stability_local, use sqrt(c/C).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_sublevel
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 11/10/2026: Initial coding;

t_start = tic;

% Check the inputs
if nargin<4 || isempty(h)
    h = {};
end
if nargin<5 || isempty(settings)
    settings = PIESOS_settings();
end
if isa(PIE,'pde_struct')
    PIE = convert(PIE,'silent');
end
if isa(g,'polyopvar')
    g = {g};
end
if isa(h,'polyopvar')
    h = {h};
end
if ~isa(PIE,'npie_struct')
    error("The PIE must be specified as 'npie_struct' (or 'pde_struct') object.")
elseif ~isa(V,'polyopvar')
    error("The functional V must be specified as 'polyopvar' object.")
elseif ~iscell(g) || isempty(g) || ~iscell(h) || ~all(cellfun(@(f) isa(f,'polyopvar'),[g(:);h(:)]))
    error("g and h must be specified as 'polyopvar' objects, or cells of 'polyopvar' objects (g nonempty).")
end
ng = numel(g);
nh = numel(h);

% The values at x=0 of V, g_i, h_j
x = PIE.x;
x0 = arrayfun(@(d) zeros(d,1),V.dims.in,'UniformOutput',false);
a = zeros(ng,1);
for i=1:ng
    a(i) = double(apply_poly(g{i},x0));
end
if any(a<0)
    error("The inequalities must hold at x=0: g_i(0) >= 0.")
elseif abs(double(apply_poly(V,x0)))>0 || any(cellfun(@(f) abs(double(apply_poly(f,x0)))>0,h))
    error("V and the equalities h_j must vanish at x=0.")
end

% Options of the SOS constraint, as for the derivative condition in
% PIESOS_stability_local
sos_opts = struct();
if settings.sos.psatz(3)>0
    sos_opts.psatz = settings.sos.psatz(3);
end
sos_opts.degs = max(settings.sos.mon_deg(3),settings.LF.mon_deg+2);     % NaN ignored by max
sol_opts.simplify = settings.simplify;
if ~isempty(settings.solver)
    sol_opts.solver = settings.solver;
end

% For each inequality, minimize p_i >= 0 such that (g_i-a_i) + p_i*V +
% sum_j q_ij*h_j is SOS, for free multipliers q_ij
info = struct('feasible',true,'c',Inf,'c_i',NaN(ng,1));
prog_sol = [];
for i=1:ng
    dpvar p
    q = [];
    if nh>0
        q = dpvar(arrayfun(@(j) sprintf('q%d',j),(1:nh)','UniformOutput',false));
    end
    prog = piesos_program(x,[p;q]);
    prog = piesos_setobj(prog,p);
    prog = piesos_ineq(prog,p);
    cons = (g{i} - a(i)) + p*V;
    for j=1:nh
        cons = cons + q(j)*h{j};
    end
    prog = piesos_ineq(prog,cons,sos_opts);
    prog_sol = piesos_solve(prog,sol_opts);
    sinfo = prog_sol.solinfo.info;
    if sinfo.pinf || sinfo.dinf || sinfo.numerr || abs(sinfo.feasratio-1)>settings.tol_feas
        info.feasible = false;
    else
        p_sol = double(piesos_getsol(prog_sol,p));
        info.c_i(i) = a(i)/max(p_sol,0);
        if a(i)==0                      % origin on the boundary of g_i >= 0
            info.c_i(i) = 0;
        end
    end
    if settings.verbose
        fprintf(" --- g_%d: c_%d = %g (%.0f s elapsed) ---\n",i,i,info.c_i(i),toc(t_start))
    end
end

% Extract the solution
if info.feasible
    info.c = min(info.c_i);
else
    info.c = 0;
end
info.time = toc(t_start);
if settings.verbose
    if info.feasible
        fprintf(" --- Sublevel set V <= %g verified in the set G ---\n",info.c)
    else
        fprintf(" --- No sublevel set could be verified in the set G ---\n")
    end
end

end
