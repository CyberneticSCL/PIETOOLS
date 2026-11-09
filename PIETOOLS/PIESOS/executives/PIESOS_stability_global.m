function [prog_sol,info] = PIESOS_stability_global(PIE,lambda,LF_form,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG_SOL,INFO] = PIESOS_stability_global(PIE,LAMBDA,LF_FORM,SETTINGS)
% tests global exponential stability, at decay rate LAMBDA, of the
% (nonlinear) PIE with the PDE state u=T*x,
%       d/dt (T*x(t)) = f(x(t)).
% It searches for a distributed polynomial Lyapunov functional V(x), of the
% form <Zop*Z(.),Pmat*Zop*Z(.)> with Z(.)=[.,...,.^d] and Pmat>=0, in
% the PIE states x or in the PDE states u (see LF_FORM), such that
%   V >= eppos*||u||^2,
% and such that, for all x, the derivative dV of V along solutions to the
% PIE satisfies
%   dV + 2*lambda*V <= 0.
% This is enforced as a distributed SOS constraint. Then V(t) decays as
% exp(-2*lambda*t), and so does ||u(t)||^2 <= V(t)/eppos. If V is
% quadratic (settings.LF.dist_deg=1), also an upper bound V <= C*||u||^2
% can be enforced, with the constant C minimized or specified (see
% settings.LF.C), so that
%   ||u(t)|| <= sqrt(C/eppos)*exp(-lambda*t)*||u(0)||.
% For settings.LF.lower='x_pde', the lower bound is instead
%   V >= eppos*||x_pde||^2,
% with x_pde=X_pde*T*x the actual PDE state (PIE.x_pde), which for e.g.
% delays leaves out the history states in u, whereas the upper bound
% remains in terms of u. Then, only the actual state is bounded,
%   ||x_pde(t)|| <= sqrt(C/eppos)*exp(-lambda*t)*||u(0)||.
%
% INPUT
% - PIE:        'npie_struct' object representing the PIE, as output by
%               `convert` for a nonlinear PDE. A 'pde_struct' object is
%               converted to a PIE;
% - lambda:     (optional) nonnegative scalar, the decay rate. Defaults
%               to 0, testing stability;
% - LF_form:    (optional) 'pie' or 'pde', the form of the Lyapunov
%               functional V: a functional in the PIE states x, with the
%               lower bound V >= eppos*||u||^2 enforced as an SOS
%               constraint ('pie'), or V = (SOS polynomial in the PDE
%               states u=T*x) + eppos*||u||^2 ('pde'). Defaults to 'pie';
% - settings:   (optional) 'struct' of settings, as output by
%               `PIESOS_settings` (see there for the fields used). Only
%               settings.LF, sos, and the solver settings are used.
%               Defaults to PIESOS_settings();
%
% OUTPUT
% - prog_sol:   PIESOS program structure of the solved program;
% - info:       'struct' with fields
%       feasible:   logical, true if the program was solved, so that the
%                   test is passed;
%       solinfo:    'struct' with the information of the solver;
%       V, dV:      'polyopvar' objects representing the solved Lyapunov
%                   functional and its derivative ([] if not solved);
%       C:          scalar, the constant C in the upper bound on V: the
%                   optimal value if settings.LF.C=[], settings.LF.C if it
%                   is a scalar, and Inf if no bound is enforced (NaN if
%                   not solved or LF.dist_deg>1);
%       gain:       scalar, sqrt(C/eppos), the gain of the exponential
%                   bound on ||u(t)|| (or ||x_pde(t)|| for LF.lower=
%                   'x_pde'; NaN if C is);
%       lambda:     the decay rate;
%       time:       the time taken, in seconds;
%
% NOTES
% - If settings.LF.dist_deg>1, V has terms of degree higher than 2, so that
%   there is no global upper bound V<=C*||u||^2: it is not enforced.
% - For settings.LF.C=Inf, no upper bound is enforced, and the program is a
%   feasibility problem, without objective.
% - For settings.LF.lower='x_pde', PIE.f_pde=X_pde*f must be computable.
%   It is not if an ODE state was declared after the PDE states, so that
%   X_pde permutes the states (tensopmat composition error), in which case
%   x_pde and u have the same norm, and 'u' is equivalent.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_stability_global
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
% DJ (CC), 11/09/2026: Initial coding;

t_start = tic;

% Check the inputs
if nargin<2 || isempty(lambda)
    lambda = 0;
end
if nargin<3 || isempty(LF_form)
    LF_form = 'pie';
end
if nargin<4 || isempty(settings)
    settings = PIESOS_settings();
end
if isa(PIE,'pde_struct')
    PIE = convert(PIE,'silent');
end
states = lower(char(LF_form));
LF = settings.LF;
lb_state = char(LF.lower);      % state in the lower bound: 'u' or 'x_pde'
eppos = LF.eppos;
C_set = LF.C;
if ~isa(PIE,'npie_struct')
    error("The PIE must be specified as 'npie_struct' (or 'pde_struct') object.")
elseif ~isscalar(lambda) || lambda<0
    error("The decay rate must be a nonnegative scalar.")
elseif ~ismember(states,{'pie','pde'})
    error("The form of the Lyapunov functional must be 'pie' or 'pde'.")
elseif ~ismember(lb_state,{'u','x_pde'})
    error("settings.LF.lower must be 'u' or 'x_pde'.")
elseif ~isempty(C_set) && ~(isnumeric(C_set) && isscalar(C_set) && C_set>0)
    error("settings.LF.C must be empty, Inf, or a positive scalar.")
elseif strcmp(lb_state,'x_pde') && isempty(PIE.space_pde.x.dims)
    error("settings.LF.lower='x_pde' requires a PIE constructed from a PDE.")
end
% Without a ball, there is no global upper bound on V of degree above 2
use_upper = LF.dist_deg==1 && ~isequal(C_set,Inf);
min_C = use_upper && isempty(C_set);        % minimize C, else C is fixed
if LF.dist_deg>1 && settings.verbose
    disp(" --- No global upper bound on V for LF.dist_deg>1: not enforced ---")
end
% ODE states, without spatial variables, only have a multiplier term
has_ode = any(~any(PIE.space.x.depmat,2));

% Extract the PIE state x, the PDE state u=T*x, and the state in the lower
% bound on V: u, or the actual PDE state x_pde=X_pde*T*x, which leaves out
% e.g. the history states of delays
x = PIE.x;
f = PIE.f;
u = PIE.T*x;
u_low = u;
f_low = f;
if strcmp(lb_state,'x_pde')
    u_low = PIE.x_pde;
    f_low = PIE.f_pde;
end

% Initialize the program, minimizing the constant C in the upper bound on V
% if it is not specified
dpvar C
if min_C
    prog = piesos_program(x,C);
    prog = piesos_setobj(prog,C);
else
    prog = piesos_program(x);
    C = C_set;
end

% Declare the Lyapunov functional V, and its derivative
Popts.degs = LF.mon_deg;
if ischar(LF.exclude)
    % Multiplier operators only for V in the PDE states, all terms for
    % the PIE states
    Popts.exclude = [0;1;1]*strcmp(states,'pde');
else
    Popts.exclude = LF.exclude;
end
Popts.sep = LF.sep;
if LF.psatz>0
    Popts.psatz = LF.psatz;
end
Popts.pos = strcmp(states,'pde');
[prog,V,dV] = piesos_LF(prog,PIE,states,(1:LF.dist_deg)',Popts);
if strcmp(states,'pde')
    % V = SOS polynomial in u, plus eppos*||u_low||^2, with derivative
    % 2*eppos*<u_low,f_low(x)>, so that V >= eppos*||u_low||^2
    V = V + eppos*(u_low'*u_low);
    dV = dV + 2*eppos*(u_low'*f_low);
end
if settings.verbose
    fprintf(" --- Lyapunov functional declared (%.0f s elapsed) ---\n",toc(t_start))
end

% The polynomials to be nonnegative, with the index k of the entry of the
% SOS settings (mon_deg, psatz) to use for each:
%   the lower bound V - eppos*||u_low||^2 >= 0, for V in the PIE states
%       (for V in the PDE states, it is satisfied by construction);
%   the upper bound C*||u||^2 - V >= 0, for V quadratic;
%   the decay -dV - 2*lambda*V >= 0, for all x
cons = cell(0,2);
if strcmp(states,'pie')
    cons(end+1,:) = {V-eppos*(u_low'*u_low),1};
end
if use_upper
    cons(end+1,:) = {C*(u'*u)-V,2};
end
cons(end+1,:) = {-dV-2*lambda*V,3};

% Enforce the constraints, using piesos_ineq, or manually, by declaring a
% distributed SOS polynomial W, and enforcing that it equals the polynomial
sos = settings.sos;
for c=1:size(cons,1)
    k = cons{c,2};
    sos_opts = struct();
    if sos.psatz(k)>0
        sos_opts.psatz = sos.psatz(k);
    end
    if ~isnan(sos.mon_deg(k))
        sos_opts.degs = sos.mon_deg(k);
    end
    if sos.use_ineq
        prog = piesos_ineq(prog,cons{c,1},sos_opts);
    else
        % W has a basis of monomials of half the degree of those in the
        % polynomial, in the PIE states x only
        [~,idcs_x] = ismember(x.varname,cons{c,1}.varname);
        Z_degmat = unique(floor(cons{c,1}.degmat(:,idcs_x)./2),'rows');
        Z_degmat = Z_degmat(sum(Z_degmat,2)>0,:);
        sos_opts.exclude = sos.exclude;
        sos_opts.exclude(1) = sos.exclude(1) && ~has_ode;
        [prog,W] = piesos_sosvar(prog,dmonomials(x,Z_degmat),sos_opts);
        prog = piesos_eq(prog,cons{c,1}-W);
    end
end
if settings.verbose
    fprintf(" --- SOS constraints declared (%.0f s elapsed) ---\n",toc(t_start))
end

% Solve the program
sol_opts.simplify = settings.simplify;
if ~isempty(settings.solver)
    sol_opts.solver = settings.solver;
end
prog_sol = piesos_solve(prog,sol_opts);

% Extract the solution. The feasibility ratio is only informative if C is
% minimized (an objective)
sinfo = prog_sol.solinfo.info;
feasible = ~(sinfo.pinf || sinfo.dinf || sinfo.numerr || (min_C && abs(sinfo.feasratio-1)>settings.tol_feas));
info = struct('feasible',feasible,'solinfo',sinfo,'V',[],'dV',[],'C',NaN,'gain',NaN,'lambda',lambda);
if feasible
    info.V = piesos_getsol(prog_sol,V);
    info.dV = piesos_getsol(prog_sol,dV);
    if min_C
        info.C = double(piesos_getsol(prog_sol,C));
    elseif use_upper
        info.C = C_set;
    elseif LF.dist_deg==1
        info.C = Inf;               % no upper bound
    end
    info.gain = sqrt(info.C/eppos);
end
info.time = toc(t_start);
if settings.verbose
    if feasible
        fprintf(" --- Global exponential stability verified, with decay rate %g ---\n",lambda)
    else
        fprintf(" --- Global exponential stability with decay rate %g could not be verified ---\n",lambda)
    end
end

end
