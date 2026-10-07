
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_Burgers_Local.m
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
% CR, 09/01/2026: Initial coding

% Local stability test script.

clear all; clear stateNameGenerator; close all; clc;

% 1. Define PDE and local stability parameters.

% Burgers Equation.
pvar s t
dom = [0,1];
x   = pde_var(s,dom);
v   = 1.0;
r   = pi^2-0.1;
PDE = [diff(x,t)==v*diff(x,s,2)+r*x-x*diff(x,s);
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];

n = 2; % degree of PDE.

% Treated as eppos^2, the lower bound on SOS LF.
eppos = 0.1;

% (n+1)-dim array containing parameters of weighted Sobolev ball.
alpha = [1, 0, 0]; % [1,0,0] for L2 ball.

% radius of local ball - we should have stability for any rad>0.
rad = 1.0;

% exponential decay rate.
lambda = 0;

% Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers 
% (respectively). Degree will be doubled when converted from quadratic 
% to linear format.
dist_degs = [1, 1, 1];

% Declare monomial degrees in independent variables used to respectively parametrize
% LF, p1, p2 multipliers. Optionally include mon_degs(4)=sos3_mon.
mon_degs = [3, 0, 0];

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
