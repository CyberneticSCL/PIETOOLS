
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

% radius of local ball.
r = 4.01;

% (n+1)-dim array containing parameters of weighted Sobolev ball.
alpha = [1, 0, 0];

% lower bound on SOS LF.
eppos = 1;

% exponential decay rate.
lambda = 0;

% Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers (respectively).
% Degree will be doubled when converted from quadratic to linear form.
dist_degs = [1, 2, 1];

% Declare monomial degrees in independent variables used to parametrize SOS LF and p1, p2 multipliers (respectively).
mon_degs = [4, 2, 2];


% Run local stability test.
[g,p1] = LocalStability(PDE, r, alpha, eppos, lambda, dist_degs, mon_degs);
