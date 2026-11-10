
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
% DJ, 10/28/2026: Use 'LocalStability_wave' instead of 'LocalStability_v2';

% Regional stability test script for the semilinear damped wave equation
% with Robin boundary feedback (cf. `Examples_Fridman.tex`, Sec. 4),
%   z_tt = z_ss - 2*mu*z_t + a*z + z^3,   s in [0,1],
%   z(0) = 0,   z_s(1) = -k*z(1),
% in the states z and z_t, using `LocalStability_wave`, on the ball
%   ||z_s||^2 + ||z_t||^2 <= rad^2.
% The nonlinearity f(z)=z^3 is only sector-bounded on a ball, since
% |z|<=||z_s||_{L2} for z(0)=0.
%
% The solution of the SOS program is returned as res = {C,rho,V,dV}, see
% `LocalStability_wave`: every solution with ||(z_s,z_t)(0)|| < rho
% remains in the ball and converges to zero.
%
% STATUS (10/28/2026), MOSEK, mu=2, a=5, k=3: with C fixed to 10 (MOSEK
% stalls when C is minimized), rad=0.3, eppos=0.05, lambda=0.05, the
% program (380k x 28k, 8 GB, 9 min) is solved to optimality: rho=0.0212.
% With the multiplier p*g of `LocalStability_v2` and V quadratic, it is
% infeasible (see the notes of `LocalStability_wave`).

clear all; clear stateNameGenerator; close all; clc;
% 1. Parameters of the wave equation and of the regional stability test.
mu  = 2.0;      % damping
a   = 5.0;      % open-loop unstable for k=0, a>pi^2/4
k   = 3.0;      % feedback gain, stabilizing for k>1.754 (Robin)

rad    = 0.3;   % radius of the energy ball
eppos  = 0.05;  % V >= eppos*N on the ball
lambda = 0.05;  % exponential decay rate, lambda>0 is needed for the cubic term of dV
C_fix  = 10;    % fix C to test feasibility, [] to minimize C

% 2. Run regional stability test.
res = LocalStability_wave(mu, a, k, rad, eppos, lambda, C_fix);

if ~isempty(res)
    C   = res{1}
    rho = res{2}
    V   = res{3};
    dV  = res{4};
end