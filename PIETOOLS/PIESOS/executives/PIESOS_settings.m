function settings = PIESOS_settings(preset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS = PIESOS_settings(PRESET) returns the settings of the PIESOS
% executives, PIESOS_stability_global and PIESOS_stability_local, as a
% 'struct'. PIESOS_sublevel only uses LF.mon_deg, sos.mon_deg(3), sos.psatz(3)
% and the solver fields. Change the fields of the output to adjust the
% settings, e.g.
%   settings = PIESOS_settings('light');
%   settings.LF.mon_deg = 4;
%   settings.sos.use_ineq = false;
%
% INPUT
% - preset:     (optional) char, 'light', 'default' or 'heavy', setting the
%               degrees of the monomials below. Defaults to 'default';
%
% OUTPUT
% - settings:   'struct' with fields
%
%   Lyapunov functional V (see piesos_LF), of degree 2*LF.dist_deg. The
%   form of V, in the PIE states x or in the PDE states u=T*x, is not set
%   here, but specified as input LF_form of the executives
%     LF.dist_deg   positive integer, maximal degree d of the distributed
%                   monomials Z(u)=[u,...,u^d] of V. V has degree 2*d.
%                   If d>1, there is no global upper bound on V;
%     LF.mon_deg    nonnegative integer, degree of the monomials in the
%                   spatial variables in the kernels of the basis operators
%                   Zop in V (opts.degs of piesos_LF);
%     LF.exclude    3 x 1 logical, the terms to exclude from Zop: the
%                   multiplier (1), the lower integral (2) and the upper
%                   integral (3) term, or 'auto'. For 'auto', [0;0;0] is
%                   used if LF_form is 'pie', and [0;1;1], only allowing
%                   a multiplier operator, if LF_form is 'pde';
%     LF.sep        logical, true to enforce separable kernels in Zop;
%     LF.psatz      nonnegative integer, order of the psatz multipliers in
%                   the SOS matrix defining V, as in piesos_sosvar;
%     LF.eppos      positive scalar, the weight of the lower bound
%                   V >= eppos*N(u) on V, with N the norm of the test;
%     LF.C          [], Inf or positive scalar, the constant C in the upper
%                   bound V <= C*N(u) on V. If [], the bound is enforced
%                   with C minimized, if a scalar c, with C=c, and if Inf,
%                   no bound is enforced (global test only);
%     LF.lower      'u' or 'x_pde', the state in the lower bound on V. For
%                   'u', V >= eppos*N(u) with u=T*x the full state, as in
%                   the upper bound. For 'x_pde', V >= eppos*N(x_pde) with
%                   x_pde=PIE.x_pde=X_pde*T*x the actual PDE state, which
%                   leaves out e.g. the history states of delays, whereas
%                   the upper bound remains in terms of u;
%
%   SOS multipliers p1, p2 of the local test (see piesos_sosvar), p1 in
%   the upper bound on V, p2 in the derivative condition on dV. They are
%   not used by the global test
%     mult.dist_deg 1 x 2 array of nonnegative integers, the degrees of
%                   the distributed monomial bases Z(u)=[u,...,u^d] of p1
%                   and p2 (the multipliers then have degree 2*d). If 0,
%                   the multiplier is not used. p1 is only needed if
%                   LF.dist_deg>1;
%     mult.mon_deg  1 x 2 array, degrees of the monomials in the spatial
%                   variables in the kernels of the operators of p1, p2;
%     mult.psatz    1 x 2 array, orders of the psatz multipliers in p1, p2;
%     mult.exclude  3 x 1 logical, the terms to exclude from the operators
%                   of p1, p2, as LF.exclude. By default, [1;0;0] only
%                   allows integral operators;
%
%   SOS constraints on V (lower bound, if not satisfied by construction),
%   on the upper bound C*N(u)-V-p1*g, and on the derivative -dV-2*lambda*V-
%   p2*g. For the global test, only the derivative condition and (if
%   LF.dist_deg=1) the upper bound are enforced, without p1, p2, g
%     sos.use_ineq  logical, true to enforce the constraints using
%                   piesos_ineq, and false to construct them manually,
%                   declaring a distributed SOS polynomial W with
%                   piesos_sosvar, and enforcing that W is equal to the
%                   polynomial with piesos_eq. The manual construction
%                   matches piesos_ineq for LF_form='pde', but was not
%                   feasible for LF_form='pie' (tested for the heat
%                   equation with a Sobolev norm);
%     sos.mon_deg   1 x 3 array, degrees of the monomials in the spatial
%                   variables in the kernels of the operators of the SOS
%                   polynomials of the lower bound, upper bound and
%                   derivative condition, respectively. NaN to let
%                   piesos_ineq select the degree (use_ineq only);
%     sos.psatz     1 x 3 array, orders of the psatz multipliers in the
%                   SOS polynomials, as in piesos_sosvar;
%     sos.exclude   3 x 1 logical, the terms to exclude from the operators
%                   of the SOS polynomials, as LF.exclude (only used if
%                   use_ineq is false, otherwise they follow from the
%                   kernels of the polynomial);
%
%   Solver
%     solver        char, the SDP solver (e.g. 'mosek', 'sedumi'), or ''
%                   to use the default of piesos_solve;
%     simplify      logical, true to simplify the SDP before solving;
%     tol_feas      scalar, if C is minimized, the program is taken to be
%                   solved if the feasibility ratio of the solver is within
%                   tol_feas of 1 (and the program is not infeasible). If C
%                   is not minimized, the ratio is not informative and not
%                   used;
%     verbose       logical, true to display the progress;
%
% NOTES
% - The defaults of the presets 'light', 'default' and 'heavy' differ in
%   the degrees LF.mon_deg, mult.mon_deg and sos.mon_deg only.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_settings
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

if nargin<1 || isempty(preset)
    preset = 'default';
end

% Distributed monomials and operator degrees of V
settings.LF.dist_deg = 1;
settings.LF.mon_deg = 3;
settings.LF.exclude = 'auto';
settings.LF.sep = false;
settings.LF.psatz = 0;
settings.LF.eppos = 1;
settings.LF.C = [];
settings.LF.lower = 'u';        % lower bound in u=T*x, or in x_pde

% SOS multipliers of the local test
settings.mult.dist_deg = [0,1];
settings.mult.mon_deg = [3,3];
settings.mult.psatz = [0,0];
settings.mult.exclude = [1;0;0];

% SOS constraints
settings.sos.use_ineq = true;
settings.sos.mon_deg = [3,3,4];
settings.sos.psatz = [1,1,1];
settings.sos.exclude = [1;0;0];

% Solver
settings.solver = '';
settings.simplify = true;
settings.tol_feas = 0.1;
settings.verbose = true;

% Degrees of the presets
switch lower(preset)
    case 'light'
        settings.LF.mon_deg = 2;
        settings.mult.mon_deg = [2,2];
        settings.sos.mon_deg = [2,2,3];
    case 'default'
        % As above
    case 'heavy'
        settings.LF.mon_deg = 4;
        settings.mult.mon_deg = [4,4];
        settings.sos.mon_deg = [4,4,5];
    otherwise
        error("Preset must be specified as 'light', 'default' or 'heavy'.")
end

end
