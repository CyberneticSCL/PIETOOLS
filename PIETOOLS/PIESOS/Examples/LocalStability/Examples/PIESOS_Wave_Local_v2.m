
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_Wave_Local_v2
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
% DJ, 10/21/2026: Initial coding;

% Regional stability test script for the semilinear damped wave equation
% with boundary feedback (cf. `Examples_Fridman.tex`, Sec. 4),
%   z_tt = z_ss - 2*mu*z_t + a*z + z^p,   s in [0,1]  (p=nl_pow, 3 or 2),
%   z(0) = 0,   z_s(1) = -k*z(1)            (feedback='robin'), or
%   z(0) = 0,   z_s(1) = -k*int_0^1 z ds    (feedback='nonlocal'),
% in the states u1=z, u2=z_t, using `LocalStability_v2`, on the ball
%   ball='energy':  ||z_s||^2 + ||z_t||^2 <= rad^2,
%   ball='zs':      ||z_s||^2             <= rad^2.
% The nonlinearity f(z)=z^p is only sector-bounded on a ball, since
% |z|<=||z_s||_{L2} for z(0)=0.
%

clear all; clear stateNameGenerator; close all; clc;
% 1. Define PDE and regional stability parameters.

% Semilinear damped wave equation, u1 = z, u2 = z_t.
pvar s t
dom = [0,1];
mu  = 2.0;      % damping
a   = 5.0;      % open-loop unstable for k=0, a>pi^2/4
k   = 3.0;      % feedback gain, stabilizing for k in (1.908,6.585) for a=5
                % (nonlocal) or k>1.754 (Robin)
feedback = 'robin';     % 'robin' or 'nonlocal'
nl_pow = 3;             % nonlinearity f(z)=z^nl_pow
u2_bc = false;          % impose u2(0)=0, with u2 an H^1 state (u2=pde_var(s,dom,1))
u1  = pde_var(s,dom);
if u2_bc
    u2 = pde_var(s,dom,1);
else
    u2 = pde_var(s,dom);
end
if strcmp(feedback,'robin')
    BC_right = subs(diff(u1,s),s,dom(2)) + k*subs(u1,s,dom(2)) == 0;
else
    BC_right = subs(diff(u1,s),s,dom(2)) + k*int(u1,s,dom(1),dom(2)) == 0;
end
PDE = [diff(u1,t) == u2;
       diff(u2,t) == diff(u1,s,2) - 2*mu*u2 + a*u1 + u1^nl_pow;
       subs(u1,s,dom(1)) == 0;
       BC_right];
if u2_bc
    PDE = [PDE; subs(u2,s,dom(1)) == 0];
end

n = 2; % degree of PDE.

% Weights of the weighted Sobolev norms, with rows in the order the PDE
% states were declared (u1=z, u2=z_t) and column j+1 for the jth derivative.
%   N_alpha = ||z_s||^2 + ||z_t||^2 (energy) or ||z_s||^2 (zs)
ball = 'energy';    % 'energy' or 'zs'
Enorm = [0,1,0; 1,0,0];
if strcmp(ball,'energy')
    alpha = Enorm;
    eppos = 0.05;   % V >= eppos*N_alpha, scalar --> eppos*alpha.
    M     = 1;      % V <= C*M*N_alpha, scalar --> M*alpha.
else
    % V must also be bounded above and below in terms of ||z_t||, that does
    % not appear in the ball: use the energy norm for the bounds
    alpha = [0,1,0; 0,0,0];
    eppos = 0.05*Enorm;
    M     = Enorm;
end

% Radius of the local ball.
rad = 0.3;

% Exponential decay rate. lambda>0 is needed for the cubic term of dV to be
% dominated (V strict Lyapunov).
lambda = 0.05;

% Degrees of the distributed monomial bases of V, p1 and p2. V is
% quadratic, so no multiplier p1 is needed in the global upper bound.
dist_degs = [1, 1, 1];

% Degrees of the spatial monomials in the operators of V, p1, p2 and sos3.
% mon_degs(2)=0 --> no multiplier p1 (not needed).
% mon_degs(1)=2 is the smallest degree for which the linear problem is feasible.
mon_degs = [2, 0, 1, 2];

% Psatz orders for sos1, sos2, sos3.
psatz = [0, 0, 0];

% Conditions to test - all must be true for full stability test.
active_conds = [true, true, true];

% Degrees of the monomials in the operators of sos1, sos2, sos3 (NaN -->
% selected from the kernels, or mon_degs(4) for sos3).
sos_degs = [NaN, NaN, NaN];


% 2. Run regional stability test.
res = LocalStability_v2(PDE, rad, alpha, eppos, lambda, dist_degs, ...
                        mon_degs, psatz, active_conds, sos_degs, M);

if ~isempty(res)
    C     = res{1}
    M_sol = res{2}
    V     = res{3};
    dV    = res{4};
    rho   = res{5}
end
