function [prog_sol,info] = PIESOS_stability_local(PIE,r,alpha,lambda,LF_form,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG_SOL,INFO] = PIESOS_stability_local(PIE,R,ALPHA,LAMBDA,LF_FORM,
%                                                               SETTINGS)
% tests local exponential stability, at decay rate LAMBDA, of the
% (nonlinear) PIE with the PDE state u=T*x,
%       d/dt (T*x(t)) = f(x(t)),
% on the ball N(u) <= R^2 in the weighted Sobolev norm with weights ALPHA,
%   N(u) = sum_i sum_j alpha(i,j+1)*||d^j/ds^j u_i||_{L2}^2.
% It searches for a distributed polynomial Lyapunov functional V, and
% SOS multipliers p1, p2, such that, with g = R^2 - N(u) >= 0 on the ball,
%   V - eppos*N >= 0,                   (lower bound)
%   C*N - V - p1*g >= 0,                (upper bound)
%   -dV - 2*lambda*V - p2*g >= 0,       (decay on the ball)
% for the derivative dV of V along solutions to the PIE, with the constant
% C minimized, or specified by settings.LF.C. The same norm N is used in
% the lower and upper bound. Then
% every solution with N(u(0)) < rho^2, with rho = R*sqrt(eppos/C), remains
% in the ball and satisfies
%   ||u(t)||_N <= sqrt(C/eppos)*exp(-lambda*t)*||u(0)||_N.
% Here, V(u(0)) <= C*N(u(0)) < eppos*R^2, which is the minimal value of V
% on the boundary N=R^2 of the ball, so that the sublevel set of V
% containing u(0) lies in the ball, and V decays along solutions.
%
% INPUT
% - PIE:        'npie_struct' object representing the PIE, as output by
%               `convert` for a nonlinear PDE. A 'pde_struct' object is
%               converted to a PIE. Coupled systems, with several PDE
%               states (possibly vector-valued, or ODE states), are
%               supported;
% - r:          positive scalar, the radius of the ball;
% - alpha:      the weights of the norm N: either an array, specifying the
%               weights of all PDE states, or a cell with the array of
%               weights of each PDE state, in the order in which the
%               states were declared. For a state depending on N spatial
%               variables s_1,...,s_N, the array alpha(k_1+1,...,k_N+1)
%               is the weight of the squared L2 norm of the derivative of
%               order k_j in s_j, j=1,...,N (see weighted_sobolev_norm): a
%               scalar for an ODE state (N=0), a vector for N=1 (use alpha
%               = [w0,w1,...] for the weights of the derivatives of order
%               0,1,...), and an array with N dimensions for N>1.
%               Derivatives with zero weight are never taken. The same
%               weights are used for each component of a vector-valued
%               state. A component of the PIE that was obtained from a PDE
%               state but depends on other variables, e.g. the history of
%               a delayed state, gets only the L2 weight. For example,
%               alpha=1 gives the L2 norm, and alpha=[1,1] the H1 norm of
%               states depending on one variable. Defaults to 1;
% - lambda:     (optional) nonnegative scalar, the decay rate. Defaults
%               to 0, testing stability;
% - LF_form:    (optional) 'pie' or 'pde', the form of the Lyapunov
%               functional V: a functional in the PIE states x, with the
%               lower bound V >= eppos*N enforced as an SOS constraint
%               ('pie'), or a SOS polynomial in the PDE states u=T*x, plus
%               eppos*N ('pde'). The form 'pde' can only be used if the
%               weights alpha involve no derivatives. Defaults to 'pie';
% - settings:   (optional) 'struct' of settings, as output by
%               `PIESOS_settings` (see there for the fields). Defaults to
%               PIESOS_settings();
%
% OUTPUT
% - prog_sol:   PIESOS program structure of the solved program;
% - info:       'struct' with fields
%       feasible:   logical, true if the program was solved, so that the
%                   test is passed;
%       solinfo:    'struct' with the information of the solver;
%       V, dV:      'polyopvar' objects representing the solved Lyapunov
%                   functional and its derivative ([] if not solved);
%       C:          scalar, the constant C: the optimal value if
%                   settings.LF.C=[], and settings.LF.C otherwise (NaN if
%                   not solved);
%       gain:       scalar, sqrt(C/eppos), the gain of the exponential
%                   bound in the norm N (NaN if not solved);
%       rho:        scalar, the radius of the region of attraction
%                   certified in the norm N, R*sqrt(eppos/C) (NaN if not
%                   solved);
%       lambda:     the decay rate;
%       time:       the time taken, in seconds;
%
% NOTES
% - If the weights alpha involve derivatives, the lower bound V >= eppos*N
%   cannot hold for V in the PDE states only, so that LF_form='pie' must
%   be used. For LF_form='pde', V satisfies the lower bound by
%   construction.
% - The norm N is declared using `weighted_sobolev_norm`.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_stability_local
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
% DJ, CR, 11/09/2026: Initial coding, based on LocalStability_v2 by CR;

t_start = tic;

% Check the inputs
if nargin<3 || isempty(alpha)
    alpha = 1;
end
if nargin<4 || isempty(lambda)
    lambda = 0;
end
if nargin<5 || isempty(LF_form)
    LF_form = 'pie';
end
if nargin<6 || isempty(settings)
    settings = PIESOS_settings();
end
if isa(PIE,'pde_struct')
    PIE = convert(PIE,'silent');
end
states = lower(char(LF_form));
LF = settings.LF;
eppos = LF.eppos;
C_set = LF.C;
if ~isa(PIE,'npie_struct')
    error("The PIE must be specified as 'npie_struct' (or 'pde_struct') object.")
elseif ~isscalar(r) || r<=0
    error("The radius of the ball must be a positive scalar.")
elseif ~isscalar(lambda) || lambda<0
    error("The decay rate must be a nonnegative scalar.")
elseif ~ismember(states,{'pie','pde'})
    error("The form of the Lyapunov functional must be 'pie' or 'pde'.")
elseif ~isempty(C_set) && ~(isnumeric(C_set) && isscalar(C_set) && C_set>0 && isfinite(C_set))
    error("settings.LF.C must be empty or a positive scalar, not Inf: without an upper bound on V, no region of attraction is certified.")
elseif ~strcmp(LF.lower,'u')
    error("settings.LF.lower must be 'u' for the local test.")
end

% The PDE states, in the order declared, with the PDE state idx(i) from
% which the i-th component of PIE.x (block of rows of T) was obtained (0 if
% none). PIE.x is ordered differently if ODE states are declared after PDE
% states
sp = PIE.space.x;
sp_pde = PIE.space_pde.x;
idx = PIE.pde2pie.x.idx(:,1);
if isempty(idx)                 % PIE not constructed from a PDE
    idx = (1:numel(sp.dims))';
    sp_pde = sp;
end
nd = numel(sp_pde.dims);        % number of PDE states
if ~iscell(alpha)
    alpha = repmat({alpha},nd,1);       % the same weights for all states
elseif numel(alpha)~=nd
    error("alpha must specify the weights of each of the %d PDE states.",nd)
end
if any(cellfun(@(A) ~isnumeric(A) || isempty(A) || any(A(:)<0),alpha)) || ~any(cellfun(@(A) any(A(:)>0),alpha))
    error("The weights of the norm must be nonnegative, and not all zero.")
end
if any(sp.dims(idx>0)~=sp_pde.dims(idx(idx>0)))
    error("The sizes of the PIE states do not match those of the PDE states.")
end
% Weights of each component of u=T*x: those of the PDE state it was
% obtained from, if it depends on the same number of spatial variables.
% Otherwise, e.g. for the history of a delayed state, only the L2 weight is
% used, and components without PDE state get weight 0
nvars = sum(sp.depmat,2);
nvars_pde = sum(sp_pde.depmat,2);
alpha_T = repmat({0},numel(sp.dims),1);
for i=find(idx(:)'>0)
    alpha_T{i} = alpha{idx(i)};
    if nvars(i)~=nvars_pde(idx(i))
        alpha_T{i} = alpha_T{i}(1);
    end
end
% V in the PDE states cannot bound a norm with derivatives from below
if strcmp(states,'pde') && any(cellfun(@(A) any(A(2:end)),alpha_T))
    error("The weights involve derivatives: V cannot be in the PDE states only (LF_form='pde').")
end
if LF.dist_deg>1 && settings.mult.dist_deg(1)==0
    warning("V has degree above 2, but p1 is not used: the upper bound on V will be infeasible.")
end
if strcmp(states,'pie') && ~settings.sos.use_ineq
    % Manual W (Z of degree floor(deg/2), fixed exclude) was infeasible for
    % the heat equation with the weighted H1 norm, where piesos_ineq was
    % feasible (11/2026)
    warning("The manual SOS construction (sos.use_ineq=false) is unreliable for V in the PIE states: use sos.use_ineq=true.")
end
% ODE states, without spatial variables, only have a multiplier term
has_ode = any(~any(sp.depmat,2));

% Extract the PIE state x and the PDE state u=T*x, declare the norm N, and
% the ball g = r^2-u_nrm >= 0
x = PIE.x;
f = PIE.f;
u = PIE.T*x;
u_nrm = weighted_sobolev_norm(u,alpha_T);
g = r^2 - u_nrm;

% Initialize the program, minimizing the constant C in the upper bound on V
% if it is not specified
dpvar C
if isempty(C_set)
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
    % V = SOS polynomial in u, plus eppos*u_nrm, so that V >= eppos*u_nrm.
    % Without derivatives in u_nrm, its derivative is 2*<u,W*f>, with W the
    % diagonal operator of the L2 weights
    V = V + eppos*u_nrm;
    dV = dV + 2*eppos*innerprod(u,diag(repelem(cellfun(@(A) A(1),alpha_T),sp.dims(:))),f);
end
if settings.verbose
    fprintf(" --- Lyapunov functional declared (%.0f s elapsed) ---\n",toc(t_start))
end

% The polynomials to be nonnegative, with the index k of the entry of the
% SOS settings (mon_deg, psatz) to use for each. For V in the PIE states,
% the lower bound V - eppos*N >= 0 is not satisfied by construction
cons = cell(0,2);
if strcmp(states,'pie')
    cons(end+1,:) = {V-eppos*u_nrm,1};
end

% The upper bound C*u_nrm - V - p1*g >= 0 and the decay
% -dV - 2*lambda*V - p2*g >= 0 on the ball g>=0, for SOS multipliers p1, p2
% (if used), with the distributed monomial basis Z(x)=[x,...,x^d] of degree
% d=mult.dist_deg
cons(end+1,:) = {C*u_nrm-V,2};
cons(end+1,:) = {-dV-2*lambda*V,3};
mult = settings.mult;
for k=1:2
    if mult.dist_deg(k)>0
        mult_opts = struct();
        mult_opts.degs = mult.mon_deg(k);
        mult_opts.exclude = mult.exclude;
        mult_opts.exclude(1) = mult.exclude(1) && ~has_ode;
        if mult.psatz(k)>0
            mult_opts.psatz = mult.psatz(k);
        end
        [prog,p] = piesos_sosvar(prog,dmonomials(x,(1:mult.dist_deg(k))'),mult_opts);
        cons{end-2+k,1} = cons{end-2+k,1} - p*g;
    end
end

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
feasible = ~(sinfo.pinf || sinfo.dinf || sinfo.numerr || (isempty(C_set) && abs(sinfo.feasratio-1)>settings.tol_feas));
info = struct('feasible',feasible,'solinfo',sinfo,'V',[],'dV',[],'C',NaN, ...
                'gain',NaN,'rho',NaN,'lambda',lambda);
if feasible
    info.V = piesos_getsol(prog_sol,V);
    info.dV = piesos_getsol(prog_sol,dV);
    if isempty(C_set)
        info.C = double(piesos_getsol(prog_sol,C));
    else
        info.C = C_set;
    end
    info.gain = sqrt(info.C/eppos);
    info.rho = r*sqrt(eppos/info.C);
end
info.time = toc(t_start);
if settings.verbose
    if feasible
        fprintf(" --- Local exponential stability verified, with decay rate %g ---\n",lambda)
        fprintf(" --- Certified radius of attraction %g, gain %g ---\n",info.rho,info.gain)
    else
        fprintf(" --- Local exponential stability with decay rate %g could not be verified ---\n",lambda)
    end
end

end
