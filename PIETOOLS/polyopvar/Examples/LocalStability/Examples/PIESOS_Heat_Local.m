
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_Heat_Local.m
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

% Local stability test script.

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

% Treated as eppos^2, the lower bound on SOS LF.
eppos = 1.0;

% (n+1)-dim array containing parameters of weighted Sobolev ball.
alpha = [0, 1, 0]; % [0,1,0] for L2 ball on x_s.

% Radius of local ball experiments with eppos=1.0; alpha = [0, 1, 0]; lambda = 0.
% rad=0.1    with dist_degs=[1,0,1], mon_degs=[4,0,4] --> numerr= and feasratio=.
% rad=0.1   failed with dist_degs=[1,0,1], mon_degs=[3,0,4] --> numerr=2 and feasratio=-0.11.
% rad=0.1   failed with dist_degs=[1,0,1], mon_degs=[3,0,3] --> numerr=2 and feasratio=0.02.
% rad=0.1   failed with dist_degs=[1,0,1], mon_degs=[3,0,0] --> numerr=2 and feasratio=-0.36.
rad = 0.1;

% exponential decay rate.
lambda = 0;

% Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers 
% (respectively). Degree will be doubled when converted from quadratic 
% to linear format.
dist_degs = [1, 0, 1];

% Declare monomial degrees in independent variables used to respectively parametrize
% LF, p1, p2 multipliers. Optionally include mon_degs(4)=sos3_mon.
mon_degs = [4, 0, 4];

% Run local stability test.
% C can be passed in as optional final argument if it is fixed.
res = LocalStability(PDE, rad, alpha, eppos, lambda, dist_degs, mon_degs);

if ~isempty(res)
    C     = res{1}
    M     = res{2}
    V     = res{3};
    dV    = res{4};
    Pcell = res{5};
end