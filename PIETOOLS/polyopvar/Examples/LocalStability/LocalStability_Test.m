
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - LocalStability_Test.m
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

% Fisher Equation.
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
alp =  5;
bet = -1;
PDE = [diff(u,t)==diff(u,s,2) + alp*u - bet*u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];

n = 2; % degree of PDE.

% Treated as eppos^2, the lower bound on SOS LF.
eppos = 0.1;

% (n+1)-dim array containing parameters of weighted Sobolev ball.
alpha = [1, 0, 0];

% radius of local ball.
r = 4.0;

% exponential decay rate.
lambda = 0;

% Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers 
% (respectively). Degree will be doubled when converted from quadratic 
% to linear format.
dist_degs = [1, 1, 1];

% Declare monomial degrees in independent variables used to parametrize
%  SOS LF and p1, p2 multipliers (respectively).
mon_degs = [3, 3, 3];


% Run local stability test. If successful res = [C, M].
% C can be passed in as optional final argument if it is fixed.
res = LocalStability(PDE, r, alpha, eppos, lambda, dist_degs, mon_degs);
