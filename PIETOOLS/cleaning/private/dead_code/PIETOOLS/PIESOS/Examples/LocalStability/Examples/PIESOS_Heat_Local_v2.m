
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_Heat_Local_v2
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
% CR, 10/08/2026: Initial coding.
% DJ, 10/14/2026: Version for `LocalStability_v2`
% DJ, 10/15/2026: Updated to use new polyopvar functions;

% Local stability test script for the semilinear heat equation
%   x_t = a0*x_ss + x^2 - (k/L)*int_0^L x ds,   x_s(0) = 0,   x(L) = 0,
% on s in [0,L], with L = 1, a0 = 1 and k = 2*a0*pi^2/L^2, using
% `LocalStability_v2`. The PDE, the radius rad of the local ball, the
% weights alpha of the weighted Sobolev norm defining the ball, the
% weights eppos and M of the lower and upper bound on the Lyapunov
% functional, the decay rate lambda, the degrees dist_degs, mon_degs and
% sos_degs, the psatz orders psatz, and the conditions active_conds to
% enforce can be changed. The solution is returned as
% res = {C,M_sol,V,dV}, see `LocalStability_v2`.

clear all; clear stateNameGenerator; close all; clc;

% 1. Define PDE and local stability parameters.

% Semilinear Heat Equation.
pvar s t
dom  = [0,1];
L    = dom(2);
x    = pde_var(s,dom);
a0   = 1.0;
k    = 2*a0*pi^2 / L^2;
PDE  = [diff(x,t) == a0*diff(x,s,2) + x^2 - (k/L)*int(x,s,dom(1),dom(2));
       subs(diff(x,s,1),s,dom(1))==0; subs(x,s,dom(2))==0];

n = 2; % degree of PDE.

% Lower bound on the SOS LF, V(x) >= sum_j eppos(j)*||d^j/ds^j T*x||^2: an array like alpha,
% or a scalar, for which eppos*alpha is used.
eppos = [0,1,0];

% Upper bound on the SOS LF, V(x) <= C*sum_j M(j)*||d^j/ds^j T*x||^2, with C minimised: an array
% like alpha, or a scalar, for which M*alpha is used (default 1).
M = [0,1,0];

% (n+1)-dim array containing parameters of weighted Sobolev norm for the local ball.
alpha = [0, 1, 0]; % [0,1,0] for L2 ball on x_s.

% Radius of local ball - for a0=1, k=2*a0*pi^2/L^2, should be able to prove up to rad =< \pi^2/2.
rad = 1;

% exponential decay rate.
lambda = 0;

% Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers 
% (respectively). Degree will be doubled when converted from quadratic 
% to linear format.
dist_degs = [1, 1, 1];

% Declare monomial degrees in independent variables used to respectively parametrize
% LF, p1, p2 multipliers. Optionally include mon_degs(4)=sos3_mon.
% mon_degs(4) is the degree of the monomials in the operators of sos3 in
% each factor of a distributed monomial.
mon_degs = [2, 2, 2];

% Declare psatz orders for sos1, sos2, sos3 respectively.
% psatz(i)=0 --> no psatz multipliers.
% psatz(i)=q>0 --> multipliers (s-a)^k*(b-s)^(q-k), k=0,...,q.
psatz = [0, 0, 1];

% Conditions to test - all must be true for full stability test.
active_conds = [true, true, true];

% Degrees of the monomials in the operators of sos1, sos2, sos3 respectively.
% NaN --> selected from the kernels of the functional (sos3: mon_degs(4)).
sos_degs = [NaN, NaN, NaN];

% Run local stability test.
res = LocalStability_v2(PDE, rad, alpha, eppos, lambda, dist_degs, ...
                        mon_degs, psatz, active_conds, sos_degs, M);

if ~isempty(res)
    C     = res{1}
    M_sol = res{2}
    V     = res{3};
    dV    = res{4};
end