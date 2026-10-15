
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_Fishers_Local
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
% CR, 10/07/2026: Initial coding
% DJ, 10/15/2026: Updated to use new polyopvar functions;

% Local stability test script for the Fisher equation
%   u_t = u_ss + alp*u - bet*u^2,   u(0) = u(1) = 0,   s in [0,1],
% using `LocalStability`. The PDE, the radius rad of the local ball, the
% weights alpha of the weighted Sobolev norm defining the ball, the lower
% bound eppos on the Lyapunov functional, the decay rate lambda, the
% degrees dist_degs and mon_degs, the psatz indicators psatz, and the
% conditions active_conds to enforce can be changed. The solution is
% returned as res = {C,M,V,dV}, see `LocalStability`.

clear all; clear stateNameGenerator; close all; clc;

% 1. Define PDE and local stability parameters.

% Fisher Equation.
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
alp = 5;
bet = -1;
PDE = [diff(u,t)==diff(u,s,2) + alp*u - bet*u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];

n = 2; % degree of PDE.

% Treated as eppos^2, the lower bound on SOS LF.
eppos = 1.0;

% (n+1)-dim array containing parameters of weighted Sobolev ball.
alpha = [1, 0, 0];

% Radius of local ball experiments with eppos=1.0; alpha = [1, 0, 0]; lambda = 0.
% rad=4.0   passed with dist_degs=[1,0,1], mon_degs=[3,0,4] --> numerr=0 and feasratio=1.00.
% rad=4.125 passed with dist_degs=[1,0,1], mon_degs=[4,0,5] --> numerr=0 and feasratio=1.00. ---> should fail.
% rad=4.125 failed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=2 and feasratio=0.9863.
% rad=4.25  failed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=2 and feasratio=0.9999.
% rad=4.5   failed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=2 and feasratio=0.9966.
% rad=4.1   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00. ---> should fail.
% rad=4.0   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00. ---> should pass.
% rad=3.9   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00.
% rad=3.5   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00.
% rad=2.5   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00.
% rad=1.0   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00.
% rad=0.5   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00.
% rad=0.1   passed with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr=0 and feasratio=1.00.
% rad=0.1   failed with dist_degs=[1,0,1], mon_degs=[3,0,4] --> numerr=0 and feasratio=-0.02.
% rad=0.1   failed with dist_degs=[1,0,1], mon_degs=[3,0,3] --> numerr=2.
% rad=0.1   failed with dist_degs=[1,0,1], mon_degs=[3,0,0] --> numerr=2.
rad = 4; 

% exponential decay rate.
lambda = 0;

% Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers 
% (respectively). Degree will be doubled when converted from quadratic 
% to linear format.
dist_degs = [1, 0, 1];

% Declare monomial degrees in independent variables used to respectively parametrize
% LF, p1, p2 multipliers. Optionally include mon_degs(4)=sos3_mon.
% Verified in PIESOS_2D (LocalStability + SOS_DP) with rad=4, eppos=1, alpha=[1,0,0],
% lambda=0, dist_degs=[1,0,1], psatz=[1,1,1]:
%   mon_degs=[3,0,2,6] solved (numerr=0, feasratio=1.00, C=1.0015).
%   mon_degs=[4,0,4,6] and [3,0,4,6] also solve; 
mon_degs = [3, 0, 2, 6];

% Declare psatz indicators for sos1, sos2, sos3 respectively.
% psatz(i)=0 --> Pcell=(s-a)*(b-s)*Pcell1.
% psatz(i)=1 --> Pcell=Pcell0 + (s-a)*(b-s)*Pcell1.
% psatz(i)=2 --> Pcell=Pcell.
psatz = [1, 1, 1];

% Conditions to test - all must be true for full stability test.
active_conds = [true, true, true];

% Run local stability test.
% C can be passed in as optional final argument if it is fixed.
res = LocalStability(PDE, rad, alpha, eppos, lambda, dist_degs, mon_degs, ...
                        psatz, active_conds);

if ~isempty(res)
    C     = res{1}
    M     = res{2}
    V     = res{3};
    dV    = res{4};
end

