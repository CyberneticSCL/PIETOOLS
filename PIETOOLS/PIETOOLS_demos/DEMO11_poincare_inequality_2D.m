% DEMO11_poincare_inequality_2D.m
% 2-D analogue of DEMO3_poincare_inequality.m.
%
% The Poincare constant is the smallest c such that
%   ||x|| <= c ( ||x_{s1}||^2 + ||x_{s2}||^2 )^{1/2}
% for all x(s1,s2) on [a,b]x[c,d] that are periodic in s2 and vanish at
% s1 = a,b. We pose this as an optimization program
%   min_{gam}   gam,
%      s.t.     <x,x> - gam (<x_{s1},x_{s1}> + <x_{s2},x_{s2}>) <= 0,
% where now c = sqrt(gam). As in DEMO3, we suppose x is twice differentiable
% in each variable, declare the boundary conditions as a BVP, and use
% CONVERT to find the operator H0op such that x(t) = H0op*xf(t), for a fundamental
% state xf; then H1op = d/ds1 H0op, H2op = d/ds2 H0op give x_{s1}, x_{s2}.
% The Dirichlet conditions at s1 = a,b can be solved out, but the periodic
% conditions in s2 cannot (the BVP is not uniquely solvable). They remain a
% constraint K*xf(t) = 0 on the fundamental state, which we impose using
% Finsler's lemma with a decision-variable slack X:
%   min gam  s.t.  -( H0'H0 - gam (H1'H1 + H2'H2) + X K + K' X' ) >= 0.
% The exact value is (b-a)/pi.

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Declare the BVP
% %   x(a,s2) = x(b,s2) = 0,   x(s1,c) = x(s1,d),   x_{s2}(s1,c) = x_{s2}(s1,d)
a = -1;  b = 1;  c = 0;  d = pi;
dom = [a b; c d];
pvar s1 s2
x = pde_var([s1;s2],dom,[2;2]);     % x in W^2 in both variables
BVP = [subs(x,s1,a)==0;
       subs(x,s1,b)==0;
       subs(x,s2,c)==subs(x,s2,d);
       subs(diff(x,s2),s2,c)==subs(diff(x,s2),s2,d)];

% % Convert the BVP:  x(t) = H0op*xf(t),  K*xf(t) = 0
BVP_sol = convert(BVP);
H0op = BVP_sol.T;                   % (H0op*xf(t)) = x(t)
H1op = diff(H0op,s1);               % (H1op*xf(t)) = x_{s1}(t)
H2op = diff(H0op,s2);               % (H2op*xf(t)) = x_{s2}(t)
Kop = BVP_sol.K;                    % constraint K*xf(t) = 0


% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(H0op.vars,H0op.dom);

% % Declare decision variables:
% %   gam \in \R,   X: constraint space --> space of xf
[prob,gam] = lpidecvar(prob,'gam');     % scalar decision variable
[spH,dmH] = copvar_space_list(H0op,'in');       % spaces and dimensions of xf (columns of H0op)
[spK,dmK] = copvar_space_list(Kop,'out');       % spaces and dimensions of the constraint K*xf
dX = struct('out',dmH,'in',dmK);
sX = struct('out',{spH},'in',{spK});
[prob,X] = lpivar_cdopvar(prob,dX,sX,dom,1);

% % Set inequality constraint:
% %   gam*(H1'H1 + H2'H2) - H0'H0 - X*K - K'*X' >= 0
prob = lpi_ineq_sop(prob,gam*(H1op'*H1op+H2op'*H2op)-H0op'*H0op-X*Kop-Kop'*X');

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
poincare_constant = sqrt(double(lpigetsol_sop(prob,gam)));

echo off

fprintf(['\n If successful, ',num2str(poincare_constant,4),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant on s1 in [',num2str(a),',',num2str(b),'] is known to be (b-a)/pi=',num2str((b-a)/pi,4),'.\n']);
