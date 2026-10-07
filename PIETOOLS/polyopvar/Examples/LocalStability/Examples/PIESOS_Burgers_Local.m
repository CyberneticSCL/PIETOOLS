%
% Burgers eq.
%
% PDE: u_t = u_ss + r*u - u*u_s;
% BCs: u(t,0) = u(t,1) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1];
%
% maps from fundamental state to PDE states
% u   = (Tv) = \int_0^s T_1(s,t)v(t)dt + \int_s^1 T_2(s,t)v(t)dt;
% u_s = (Rv) = \int_0^s R_1(s,t)v(t)dt + \int_s^1 R_2(s,t)v(t)dt;
% T_1(s,t)=(s-1)*t; T_2(s,t) = s*(t-1);
% R_1(s,t)=t;       R_2(s,t) = t-1;
%
% PIE: (Tv_t) = f(v) = v + (r*Tv) - (Tv)(Rv);
%
% PIE as a polyopvar: Z_1(v) = v; C = C_1 = [[C_111, C_112] 
%                                            [C_121, C_122]
%                                            [C_131, C_132]]
%
% C_111 = 1; C_112 = 1; C_121 = 1; C_122 = r*T; C_131 = -T; C_132 = R
% where nz = 1; m_1 = 3; d_1 = d_11 = 2; size(C_1) = (3,2)


clear all; clear stateNameGenerator; close all; clc;


% 1. Define PDE and local stability parameters.

% Burgers Equation.
pvar s t
dom = [0,1];
x   = pde_var(s,dom);
r   = pi^2-0.1;
PDE = [diff(x,t)==diff(x,s,2)+r*x-x*diff(x,s);
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];

n = 2; % degree of PDE.

% Treated as eppos^2, the lower bound on SOS LF.
eppos = 0.1;

% (n+1)-dim array containing parameters of weighted Sobolev ball.
alpha = [1, 0, 0];      % <-- [1,0,0] for L2 ball

% radius of local ball.
rad = 1.0;                % for Burgers', we should have stability for any r>0

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
