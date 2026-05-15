%
% Burgers-Fishers eq.
%
% PDE: x_t = a*x_ss + x^2 - k/L * \int_0^L x ds;
% BCs: x_s(t,0) = x(t,L) = 0;
%
% fundamental state: v = x_ss \in L_2[0,1] (assuming L=1);
%
% PIE: (Tv_t) = f(v) = (a - k/L * vec(T))*v + Tv*Tv;
%
% CURRENT STATE OF SCRIPT JUST CONTAINS THE REPRESENTATION OF THE PIE.
% CAN THEN BASE REST OF SCRIPT ON SOBOLEV_FISHER ONCE BUGS HAVE BEEN FIXED.

clear;  clear stateNameGenerator

% Declare the nonlinear PDE
pvar   s t s_dum
L = 1;
a = 1;
k = 1;
dom = [0,L];
x = pde_var(s,dom);

% FULL PDE
% PDE = [diff(x,t) == a*diff(x,s,2) + x^2 - (k/L)*int(x,s,dom(1),dom(2));
%        subs(diff(x,s,1),s,dom(1))==0; subs(x,s,dom(2))==0];

% TEMP PDE FOR COMPUTING T.
PDE = [diff(x,t) == a*diff(x,s,2) + x^2;
        subs(diff(x,s,1),s,dom(1))==0; subs(x,s,dom(2))==0];

% Script parameters



%% 1. Modelling PIE.

% Convert to a PIE
PIE = convert(PDE); % CANT HANDLE THE INTEGRAL TERM BUT CAN USE IT TO COMPUTE T.
Top = PIE.T;
Top_opvar = ndopvar2dopvar(Top);
f = PIE.f;
x = polyopvar(f.varname,s,dom);

% just need to add integral term to f?
% d0 = dmonomials(x,0);
% innerprod(d0,x,Top)
